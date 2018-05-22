require(MASS)
require(pomp)

## -- Utility functions ------------------------------------------------------------
anti_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

logit <- function(x){
  return(log(x/(1-x)))
}

log_titer_trans <- function(x){
    return(log2(x/10) + 2)
}

standardize <- function(vec){
  return((vec-mean(vec))/(sd(vec)))
}

dot_product <- function(x,y){
  return(sum(x*y))
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

AICc = function(k,n,LL){
  return(2*k - 2*LL + (2*k^2 + 2*k)/(n - k - 1))
}

## data generation functions ----------------------------------------------------------------------------------
make_pomp_panel <- function(i, 
                            ind_data_list, 
                            test_params = shared_params, 
                            imprinting_subtype,
                            timestep = 3, 
                            log_transform_titer = TRUE,
                            L_tol = .00005,
                            n_years_prior = 7,
                            this_subtype = "H3N2"){
  
  cat("i is ",i, "\n")
  ind_data <- ind_data_list[[i]]$y %>% arrange(obs_time)
  n_strains <- ind_data_list[[i]]$n_strains
  vis_dates <- ind_data_list[[i]]$dates
  
  if(log_transform_titer){
    ind_data <- ind_data %>% mutate_at(.vars = vars(contains("h_obs")),.funs = funs(log_titer_trans))
  }
  
  h_init <- as.numeric(ind_data %>% filter(obs_time == 0) %>% select(contains("h_obs")))
  times_ind <- ind_data$obs_time
  n.vis <- length(times_ind)
  age_ind <- unlist(ind_data_list[[i]]$demog_init["age"])
  
  test_params["n_strains"] <- n_strains
  test_params["init_age"] <- age_ind
  test_params[paste0("h_t0_s",c(1:n_strains))] <- h_init
  test_params["p_imprinted_h3"] <- unlist(ind_data_list[[i]]$imprinting_probs["h3n2"])
  test_params["p_imprinted_group_1"] <- unlist(ind_data_list[[i]]$imprinting_probs["h1n1"]) + unlist(ind_data_list[[i]]$imprinting_probs["h2n2"])
  test_params["h_baseline_individual_observed"] <- min(ind_data$h_obs_s1)
  
  ## Covariates: Flu intensity 
  date_0 <- as.Date(vis_dates[1], orgin = '1970-1-1')
  date_end <- as.Date(vis_dates[length(vis_dates)], orgin = '1970-1-1')
  
  # Start tracking flu intensity some number of years (n_years_prior) before first visit date 
  if(age > n_years_prior)
    sim_start_date = date_0 - n_years_prior*365
  }
  if(age <= n_years_prior){
    sim_start_date = date_0 - age*365
  }
  
  L_data_sub <- L_data %>% filter(full_date >= sim_start_date & full_date <= date_end & Subtype == this_subtype & !is.na(L)) %>% 
    arrange(full_date) %>% 
    mutate(obs_time = (as.numeric(as.Date(full_date, origin = '1970-1-1')) - as.numeric(as.Date(date_0, origin = '1970-1-1'))))
  prediction_range <- seq(from = min(L_data_sub$obs_time), to = max(ind_data$obs_time), by = timestep)
  
  L_pre_sub <- L_data_sub %>% 
    filter(obs_time <= 0) %>% 
    mutate(p_L = L/sum(L))
  L_pre_sub$obs_time[nrow(L_pre_sub)] <- 0

  covartable <- data.frame(
    obs_time = prediction_range,
    L = predict(smooth.spline(x=L_data_sub$obs_time, y=L_data_sub$L),
                x=prediction_range)$y
  ) %>% mutate(L = sapply(L,max, L_tol))
  
  L_pre_sub_covar <- rep.row(data.frame(t(data.frame(c(L_pre_sub$obs_time, L_pre_sub$p_L)))), nrow(covartable))
  
  # Define the covariate table that goes into the pomp object 
  covartable <- cbind(covartable, L_pre_sub_covar) 
  names(covartable) <- c("obs_time","L",paste0("epi_times_",c(1:nrow(L_pre_sub))), paste0("p_epi_",c(1:nrow(L_pre_sub))))
  
  # Define the latent states in the model 
  statenames = c(sprintf("t_last_infection_s%d",c(1:n_strains)),
                 sprintf("t_clear_s%d",c(1:n_strains)),
                 sprintf("h_s%d",c(1:n_strains)),
                 sprintf("h_last_infection_s%d",c(1:n_strains)),
                 sprintf("I_s%d",c(1:n_strains)),
                 "h_baseline_individual",
                 "d_acute",
                 "q1",
                 "q2",
                 "q",
                 "age"
  )
  
  # Name of observation variables for the pomp model 
  obsnames = c(sprintf("h_obs_s%d",c(1:n_strains)))
  
  if("obs" %in% names(ind_data)){
    ind_data %>% select(-obs)
  }
  pomp_data <- ind_data
  log_paramnames <- names(test_params)
  
  # Load the process model file (built in "C-snippets")
  source("rprocess.R")
  pomp_object <- pomp(
    data = pomp_data,
    times ="obs_time",
    t0=0,
    params=unlist(test_params),
    covar = covartable,
    tcovar = "obs_time",
    rprocess = discrete.time.sim(step.fun=rprocess, delta.t= timestep),
    dmeasure = dmeasure,
    rmeasure = rmeasure,
    obsnames = obsnames,
    statenames = statenames,
    paramnames = log_paramnames,
    initializer = init
  ) 
  
  return(list(pomp_object = pomp_object, 
              spec_params = c(test_params["init_age"],
                              test_params[paste0("h_t0_s",c(1:n_strains))],
                              test_params["p_imprinted_group_1"],
                              test_params["p_imprinted_h3"],
                              test_params["h_baseline_individual_observed"]
              )
  )
  )
}


mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000){
  smooth_fit <- loess(lp ~ parameter,span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0,length(parameter))
  weight[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp,parameter=parameter,confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(
         parameter=parameter_grid,
         smoothed=smoothed_loglik,
         quadratic=predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
       ),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), se=sqrt(se_total_squared)
  )
}

