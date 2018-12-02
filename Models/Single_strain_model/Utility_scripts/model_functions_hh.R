### simulate HPV model 
require(MASS)
require(pomp)

## -- utility functions ------------------------------------------------------------
anti_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

logit <- function(x){
  return(log(x/(1-x)))
}

log_titer_trans <- function(x){
    return(log2(x/10) + 1)
  #return(log2(x))
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

## data generation functions ----------------------------------------------------------------------------------
make_pomp_panel_hh <- function(i, 
                               hh_data_list, 
                               test_params = shared_params, 
                                imprinting_subtype,
                                timestep = 3, 
                                MAX_MEMBERS = 5,
                                log_transform_titer = TRUE,
                                L_tol = .00005,
                                n_years_prior = 5,
                                L_data = intensity_data,
                                subtype = "pH1N1"
){
  
  cat("i is ",i, "\n")
  hh_data <- hh_data_list[[i]]$y %>% arrange(obs_time)
  n_members <- hh_data_list[[i]]$n_members
  vis_dates <- hh_data_list[[i]]$dates
  starting_dates <- hh_data_list[[i]]$starting_dates
  
  if(log_transform_titer){
    hh_data <- hh_data %>% mutate_at(.vars = vars(contains("h_obs")),.funs = funs(log_titer_trans))
  }
  
  # extract initial titer
  h_init <- as.numeric(hh_data %>% filter(obs_time == 0) %>% select(contains("h_obs")))
  h_min <- as.numeric( hh_data %>% select(contains("h_obs")) %>% apply(.,2,min, na.rm=T)) 
  times_hh <- hh_data$obs_time
  n.vis <- length(times_hh)
  
  for(j in c((n_members + 1): MAX_MEMBERS)){
    new_col <- array(NA, nrow(hh_data))
    names_old <- names(hh_data)
    hh_data <- cbind(hh_data,new_col) %>% as.data.frame()
    names(hh_data) <- c(names_old, paste0("h_obs_m",j))
  }
  
  
  test_params["n_strains"] <- n_strains
  test_params["n_members"] <- n_members
  test_params[paste0("init_age_m",c(1:n_members))] <- unlist(hh_data_list[[i]]$demog_init["age"])
  test_params[paste0("init_age_m",c((n_members + 1): MAX_MEMBERS))] <- NA
  test_params[paste0("h_t0_m",c(1:n_members))] <- h_init
  test_params[paste0("h_t0_m",c((n_members + 1): MAX_MEMBERS))] <- NA
  test_params[paste0("p_imprinted_h3_m",c(1:n_members))] <- unlist(hh_data_list[[i]]$imprinting_probs["h3n2"])
  test_params[paste0("p_imprinted_h3_m",c((n_members + 1): MAX_MEMBERS))] <- NA
  test_params[paste0("p_imprinted_group_1_m",c(1:n_members))] <- unlist(hh_data_list[[i]]$imprinting_probs["h1n1"]) + unlist(hh_data_list[[i]]$imprinting_probs["h2n2"])
  test_params[paste0("p_imprinted_group_1_m",c((n_members + 1): MAX_MEMBERS))] <- NA
  test_params[paste0("h_baseline_observed_m",c(1:n_members))] <- h_min
  test_params[paste0("h_baseline_observed_m",c((n_members + 1): MAX_MEMBERS))] <- NA
  
  ## Covariates: Flu intensity 
  date_end <- as.Date(vis_dates[length(vis_dates)], orgin = '1970-1-1')
  
  L_covar <- list()
  this_subtype = subtype
  for( i in c(1:n_members)){
    sim_start_date <- as.Date(starting_dates[which(starting_dates$member == starting_dates$member[i]),]$date, origin = "1970-1-1") - n_years_prior*365
    L_data_sub <- L_data %>% filter(full_date >= sim_start_date & full_date <= date_end & Subtype == this_subtype & !is.na(L)) %>% 
      arrange(full_date) %>% 
      mutate(obs_time = (as.numeric(as.Date(full_date, origin = '1970-1-1')) - as.numeric(as.Date(starting_dates[which(starting_dates$member == starting_dates$member[i]),]$date))))

    prediction_range <- seq(from =min(L_data_sub$obs_time), to = max(L_data_sub$obs_time), by = 2)
    
    L_pre_sub <- L_data_sub %>% 
     filter(obs_time <= 0) %>% 
     mutate(p_L = L/sum(L))
    L_pre_sub$obs_time[nrow(L_pre_sub)] <- 0
   # L_pre_sub <- data.frame(
    #  obs_time = c(prediction_range_pre,0),
    #  L = predict(smooth.spline(x=L_data_sub$obs_time, y=L_data_sub$L),
    #              x=c(prediction_range_pre,0))$y
   # ) %>% mutate(p_L = L/sum(L))
    
    covartable <- data.frame(
      obs_time = prediction_range,
      L = predict(smooth.spline(x=L_data_sub$obs_time, y=L_data_sub$L),
                  x=prediction_range)$y
    ) %>% mutate(L = sapply(L,max, L_tol))
    
    L_pre_sub_covar <- rep.row(data.frame(t(data.frame(c(L_pre_sub$obs_time, L_pre_sub$p_L)))), nrow(covartable))
    covartable <- cbind(covartable, L_pre_sub_covar) 
    names(covartable) <- c("obs_time",paste0("L_m",i),paste0("epi_times_",c(1:nrow(L_pre_sub)),"_m",i), paste0("p_epi_",c(1:nrow(L_pre_sub)),"_m",i))
    L_covar[[i]] <- as.data.frame(covartable)
  }

  covartable <- L_covar[[1]]
  if(length(L_covar) > 1){
    for( i in c(2:length(L_covar))){
      covartable <- merge(covartable,L_covar[[i]])
    }
  }  
  
  if(n_members < 4){
    for( i in c((n_members + 1):4)){
      names <- c(paste0("p_epi_1_m",i),paste0("epi_times_1_m",i))
      for(k in names){
        if(!(k %in% names(covartable))){
          covartable[,k] = NA
        }
      }
    }
  }
  cov_sub <- covartable %>% select(contains("L_"))
  covartable <- cbind(covartable %>% select(-(contains("L_"))),cov_sub)
  statenames = c(sprintf("t_last_infection_m%d",c(1:MAX_MEMBERS)),
                 sprintf("t_clear_m%d",c(1:MAX_MEMBERS)),
                 sprintf("h_m%d",c(1:MAX_MEMBERS)),
                 sprintf("h_last_infection_m%d",c(1:MAX_MEMBERS)),
                 sprintf("I_m%d",c(1:MAX_MEMBERS)),
                 sprintf("q1_m%d",c(1:MAX_MEMBERS)), 
                 sprintf("q2_m%d",c(1:MAX_MEMBERS)),
                 sprintf("q_m%d",c(1:MAX_MEMBERS)),
                 sprintf("d_acute_m%d",c(1:MAX_MEMBERS)),
                 sprintf(paste0("h_baseline_individual_m",c(1:MAX_MEMBERS))),
                 sprintf("age_m%d",c(1:MAX_MEMBERS))
  )
  obsnames = c(sprintf("h_obs_m%d",c(1:MAX_MEMBERS)))
  
  pomp_data <- hh_data
  log_paramnames <- names(test_params)
  source("../Utility_scripts/rprocess_hh.R")
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
  #test_params["logit_d2_kids"] <- -100
  #test_params["logit_d2_adults"] <- -100
  
  #s6 <- simulate(pomp_object, params = test_params)
  #df_tot5 <- as.data.frame(t(states(s6))) %>% mutate(time = s6@times)
    
  #p1 <- ggplot(df_tot5, aes(x = time/365, y = q1_m2)) + geom_line() + 
    #geom_line(aes(x = time/365, y = q1_m3, col = "red")) + xlab("Time (y)") + ylab("Susceptibility") + plot_themes
  #p2 <-  ggplot(df_tot2, aes(x = time/365, y = q1_m2)) + geom_line() + 
   # geom_line(aes(x = time/365, y = q1_m3, col = "red")) + xlab("Time (y)") + ylab("Susceptibility") + plot_themes
  
  
  
  return(list(pomp_object = pomp_object, 
              spec_params = c(test_params[paste0("init_age_m",c(1:MAX_MEMBERS))],
                              test_params[paste0("h_t0_m",c(1:MAX_MEMBERS))],
                              test_params[paste0("p_imprinted_group_1_m",c(1:MAX_MEMBERS))],
                              test_params[paste0("p_imprinted_h3_m",c(1:MAX_MEMBERS))],
                              test_params[paste0("h_baseline_observed_m",c(1:MAX_MEMBERS))],
                              test_params["n_members"]
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

# N-fold rises
err_func = function(x,thres,sig,sig2){
  if(x > thres ){
    y = floor(rnorm(1,mean = x,sd = sig))
    if(y > 9){
      y = 9
    }
    if(y < 1){
      y = 0
    }
  }
  if(x <= thres){
    sigma = sig2
    y = floor(rnorm(1,mean = x,sd = sigma))
    if(y < 1){
      y = 0
    }
  }
  return(y)
}

err_func2 <- function(x,thres,sig,sig2){
  if(x < thres){
    sigma = sig2
    p_correct = pnorm(1,x, sigma) - pnorm(-100,x,sigma)
    correct = rbinom(1,1,p_correct)
    if(correct == 1){
      x_obs = 0
    }
    if(correct == 0 ){
      value_vec = c(1:11)
      y = array(NA, (length(value_vec) - 1))
      for(i in c(1:length(y))){
        p_correct = pnorm(value_vec[i+1],x, sigma) - pnorm(value_vec[i],x,sigma)
        y[i] = p_correct
      }
      x_obs = value_vec[which(rmultinom(1,1,y) == 1)]
      if(x_obs > 9){
        x_obs = 9
      }
    }
  }
  
  if(x >= thres){
    sigma = sig
    value_vec = c(-10:11)
    y = array(NA, (length(value_vec) - 1))
    for(i in c(1:length(y))){
      #cat("i is ", i, "\n")
      p_correct = pnorm(value_vec[i+1],x, sigma) - pnorm(value_vec[i],x,sigma)
      #cat("p correct is", p_correct, "\n")
      y[i] = p_correct
    }
    x_obs = value_vec[which(rmultinom(1,1,y) == 1)]
    if(x_obs < 1){
      x_obs = 0
    }
    if(x_obs > 9){
      x_obs = 9
    }
  }
  return(x_obs)
}


