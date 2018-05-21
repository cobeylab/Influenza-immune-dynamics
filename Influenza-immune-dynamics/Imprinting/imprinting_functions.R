## -----------------------------------------------------------------------------------------------------------------------------------------------
## Calculate the probability of exposure to influenza A subtypes given data on flu intensitiy, subtype frequencies, and birthdate
## -----------------------------------------------------------------------------------------------------------------------------------------------

## Note: this code calculates an individual's imprinting probabilities for the first visit in the study. 
##    For the adults analyzed in the paper (ages 35-50), the imprinting probability is constant over the observation 
##    period. If you are calculating imprinting probabilities for individuals < 12 y old, you must re-calculate the 
##    probability at every observation.

prob_primary_flu_next_dt <- function(p_uninfected = 1, 
                                     attack_rate = 0.28/12, # Using attack rate of 0.28 reported by Gostic et al. 2016
                                     flu_data = relative_incidence_full
){
  AR = attack_rate # currently no interannual variation 
  p_naive = p_uninfected*(1-AR)
  p_first_infection = p_uninfected*AR
  return(list(p_naive = p_naive, p_first_infection = p_first_infection))
}

prob_primary_exposure<- function(id,
                                 test_year = 2009,
                                 test_month = 1,
                                 strain_vec = c("h1n1","h2n2","h3n2"),
                                 p_uninfected_init = 1,
                                 incidence_data = relative_incidence,
                                 demog_data = demography_full,
                                 sero_data = serology_full,
                                 resolution = "year",
                                 attack_rate = 0.28/12,
                                 duration_mat_Abs = duration_maternal_Abs,
                                 n_yrs =12,
                                 MONTHS_PER_YEAR = 12
){
  cat("id is ",id, "\n")
  demog <- demog_data %>% filter(memberID == id)
  sero <- sero_data %>% filter(memberID == id)
  first_vis_date <- min(sero$date)
  first_vis_month <- month(first_vis_date)
  first_vis_year <- year(first_vis_date)
  birthdate <- as.Date(demog$birthdate,origin = '1970-1-1') %m+% months(duration_mat_Abs)
  birthyear <- year(birthdate)
  birthmonth <- month(birthdate)
  p_uninfected = p_uninfected_init
  df_exposure <- expand.grid(year = c(birthyear:(birthyear + (n_yrs+1))), month = c(1:12)) %>% arrange(year)

  df_exposure$p_infected = NA
  df_exposure$p_naive = NA
  
  # Make sure that the exposure only accounts for time after birth and before first visit 
  df_exposure <- df_exposure %>% filter(!(year <=birthyear & month<=birthmonth))  %>% filter(!(year>test_year)) %>% filter(!(year == first_vis_year & month > test_month))
  if(nrow(df_exposure) > n_yrs*MONTHS_PER_YEAR){
    df_exposure <- df_exposure[1:(n_yrs*MONTHS_PER_YEAR),]
  }
  
  df_exposure[1,]$p_infected <- prob_primary_flu_next_dt(p_uninfected = p_uninfected_init)["p_first_infection"]
  df_exposure[1,]$p_naive <-  prob_primary_flu_next_dt(p_uninfected = p_uninfected_init)["p_naive"]
  for(y in c(2:nrow(df_exposure))){
    df <- prob_primary_flu_next_dt(p_uninfected = unlist(df_exposure[y-1,]$p_naive))
    df_exposure[y,]$p_infected <- unlist(df["p_first_infection"])
    df_exposure[y,]$p_naive <- unlist(df["p_naive"])
  }
  df_exposure <- as.data.frame(apply(df_exposure,2,unlist))
  
  df_exposure["h1n1"] = 0 
  df_exposure["h3n2"] = 0
  df_exposure["h2n2"] = 0
  
  for(s in c(1:length(strain_vec))){
    strain = strain_vec[s]
    for(y in c(1:nrow(df_exposure))){
      if(resolution == "month"){
        strain_data <- incidence_data %>% filter(year == df_exposure[y,]$year & month == df_exposure[y,]$month) 
      }
      if(resolution == "year"){
        strain_data <- incidence_data %>% filter(year == df_exposure[y,]$year) 
      }
      p_inf_this_strain <- as.numeric(strain_data %>% select(contains(strain)))*df_exposure[y,]$p_infected
      df_exposure[y,which(names(df_exposure)==strain)] <- p_inf_this_strain
    }
  }
  cum_probs <- apply(df_exposure[, names(df_exposure) %in% strain_vec ],2,sum)
  
  ## Scaling 
  if(nrow(df_exposure) == (n_yrs*MONTHS_PER_YEAR)){
    strain_probs <- cum_probs/sum(cum_probs)
    strain_probs["p_naive"] <- 0
  }
  
  if(nrow(df_exposure) < (n_yrs*MONTHS_PER_YEAR)){
    p_naive <- df_exposure[nrow(df_exposure),]$p_naive
    strain_probs <- (cum_probs*(1 - p_naive))/sum(cum_probs)
    strain_probs["p_naive"] <- p_naive
  }
  return(strain_probs)
}
