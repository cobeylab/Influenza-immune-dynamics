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
  if(x > 0){
    return(log2(x/10) + 2)
  }
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


## Pomp data generation functions ----------------------------------------------------------------------------------
make_pomp_panel <- function(i, 
                            test_ids,
                            data, 
                            test_params = shared_params, 
                            timestep = 2.5, 
                            log_transform_titer = TRUE
                            ){
  
  cat("i is ",i, "\n")
  ind_data <- data %>% filter(memberID == test_ids[i]) %>% 
    select(date, value, swab_date, age_at_recruitment) %>% 
    rename(h_obs = value) %>% 
    mutate(obs_time = as.numeric(date - min(date)),
           swab_time = as.numeric(swab_date - min(date)))
  
  vis_dates <- ind_data$date
  
  if(log_transform_titer){
    ind_data <- ind_data %>% mutate_at(.vars = vars(contains("h_obs")),.funs = funs(log_titer_trans))
  }

  # extract initial titer
  h_init <- as.numeric(ind_data %>% filter(obs_time == 0) %>% select(contains("h_obs")))
  times_ind <- ind_data$obs_time
  n.vis <- length(times_ind)
  

  test_params["init_age"] = ind_data$age_at_recruitment
  test_params["h_t0"] = h_init
  test_params["t_infection"] = unique(ind_data$swab_time)
 
  statenames = c("h",
                 "h_init",
                 "age"
  )
  obsnames = "h_obs"

  pomp_data <- ind_data %>% select(
    obs_time,h_obs
  )
  log_paramnames <- names(test_params)
  source("rprocess_boosting_model.R")
  pomp_object <- pomp(
    data = pomp_data,
    times ="obs_time",
    t0=0,
    params=unlist(test_params),
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
                              test_params["h_t0"],
                              test_params["t_infection"]
                              )
              )
  )
}

## Calculate 95% CIs by Monte Carlo Adjusted profiles ## ------------------------------------------------------------------
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

