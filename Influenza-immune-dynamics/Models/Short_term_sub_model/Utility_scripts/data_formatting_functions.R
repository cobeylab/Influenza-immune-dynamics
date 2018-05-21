
log_titer_trans <- function(x){
  return(log2(x/10) + 2)
}

get_age_at_recruitment <- function(memberID, demography_table = demog){
  return(unique(demography_table[demography_table$memberID == memberID,]$age_at_recruitment))
}

get_birth_date <- function(memberID, demography_table = demog){
  return(unique(demography_table[demography_table$memberID == memberID,]$birthdate))
}

get_population_flu_data_for_date_range <- function(date_range, pop_data = flu_pop ){
  start <- min(date_range) 
  end <- max(date_range) 
  pop_data$date_start <- as.Date(paste(pop_data$year, pop_data$week_start, sep = "-"), format = "%Y-%m-%d")
  pop_data$date_end <- as.Date(paste(pop_data$year, pop_data$week_end, sep = "-"), format = "%Y-%m-%d")
  pop_data_subset <- pop_data %>% filter(date_start >= start & date_end <= end)
  return(pop_data_subset)
}

get_initial_age_for_ind <- function(id, dem_data = demog, DAYSPERYEAR = 365){
  data <- dem_data %>% filter(memberID == id)
  dem_df <- data.frame(member = data$memberID,
                        age = as.numeric(as.Date(data$visitdate) - as.Date(data$birthdate))/DAYSPERYEAR,
                        gender = data$gender)
  return(dem_df)
}

get_serology_data_for_ind_and_strain <- function(strain_id, id , strain_vec = strains, serology_data = serology){
  strain_name <- strain_vec[strain_id]
  data <- serology %>% filter(memberID == id & name == strain_name) %>% as.data.frame %>% select(-c(hhid,member,visit_id,name)) %>% 
    filter(!is.na(date))
  n_strains <- length(strain_vec)
  data <- data %>% select(-memberID) %>% arrange(date)
  names(data) <- c("date",paste0("h_obs_s",strain_id))
  data$obs_time <- (as.numeric(as.Date(data$date,origin = '1970-1-1')) - as.numeric(as.Date(data$date[1],origin = '1970-1-1')))
  data <- data %>% arrange(obs_time)
  data <- melt(data, id.vars = c("obs_time","date"))
  return(data)
}

get_PCR_data_for_ind_and_strain <- function(strain_id,id , strain_vec = strains, swab_data = swabs){
  strain_name <- strain_vec[strain_id]
  data <- swab_data %>% filter(memberID == id) %>% select(-c(hhid,member)) %>% arrange(date)
  data$obs_time <- (as.numeric(as.Date(data$date,origin = '1970-1-1')) - as.numeric(as.Date(data$date[1],origin = '1970-1-1')))
  eturn(data)
}

get_population_flu_data_for_date_range <- function(date_range, pop_data = flu_pop ){
  start <- min(date_range) 
  end <- max(date_range) 
  pop_data$date_start <- as.Date(paste(pop_data$year, pop_data$week_start, sep = "-"), format = "%Y-%m-%d")
  pop_data$date_end <- as.Date(paste(pop_data$year, pop_data$week_end, sep = "-"), format = "%Y-%m-%d")
  pop_data_subset <- pop_data %>% filter(date_start >= start & date_end <= end)
  return(pop_data_subset)
}

## Make pomp data
make_pomp_data_for_ind <- function(id, 
                                   strain_vec = strains, 
                                   serology_data = serology, 
                                   imprinting_data,
                                   simulate_times = F, 
                                   demographic_data = demog,
                                   n_vis_sim = NA){
  cat("id is ", id, "\n")
  n_strains <- length(strain_vec)
  hh_data_titer <- do.call("rbind",lapply(c(1:n_strains), get_serology_data_for_ind_and_strain, id = id))
  hh <- hh_data_titer %>% 
    mutate(date_diff = c(as.numeric(diff(date)),NA)) %>% 
    filter(date_diff > 20 | is.na(date_diff)) %>% 
    mutate(obs_time = obs_time - min(obs_time))
  data_all <- hh %>% select(-date_diff)
  min_titer <- min(data_all$value)
  dates <- unique(data_all$date)
  data <- reshape(data_all %>% select(-date), timevar = "variable", idvar = "obs_time", direction = "wide" )
  names(data) <- c("obs_time", as.character(unique(data_all$variable)))
  demog_init <- get_initial_age_for_ind(id = id, dem_data =demographic_data) %>% select(-member)
  imprinting_probs <- imprinting_data %>% filter(memberID == id) %>% unlist
  if(simulate_times){
    init_titer <- hh_data_titer %>% filter(obs_time == 0)
    date_0 <- dates[1]
    t_0 <- min(data$obs_time)
    t_end <- max(data$obs_time)
    obs_time <-  seq(t_0, t_end, length.out = n_vis_sim)
    data = data.frame(obs_time = obs_time)
    for(i in c(1:n_strains)){
      init_titer_val <- init_titer %>% filter(init_titer$variable == paste0("h_obs_s",i)) %>% select(value) %>% as.numeric
      data  <- cbind(data,c(init_titer_val,rep(min_titer, n_vis_sim -1)))
    }
    names(data) <- c("obs_time",paste0("h_obs_s",c(1:n_strains)))
    dates <- as.Date(as.numeric(as.Date(date_0, origin ='1970-1-1')) + obs_time, origin = '1970-1-1')
  }
  return(list(y = data, n_strains = n_strains, demog_init = demog_init, dates = dates, imprinting_probs = imprinting_probs))
} 

get_frac_missing <- function(id,df = unvacc_HHs){
  df_sub <- df %>% filter(hhid == id) 
  total_NA = sum(is.na(select(df_sub,starts_with("vax"))))/(sum(!is.na(select(df_sub,starts_with("vax")))) + sum(is.na(select(df_sub,starts_with("vax")))))
  return(total_NA)
}

get_memberID <- function(i,table1,table2 = demography_full){
  #cat("i is ",i, "\n")
  mem <- table2[which(table2$hhid == table1[i,]$hhid & table2$member == table1[i,]$member),]$memberID
  return(mem)
}

get_study_duration <- function(id, data = serology){
  data_sub <- data %>% as.data.frame() %>% filter(memberID == id) %>% select(c(visit_id,date))  %>% distinct()
  df <- data.frame(memberID = id,
                   n_vis = max(data_sub$visit_id),
                   duration = max(as.Date(data_sub$date,origin = "1970-1-1")) - min(as.Date(data_sub$date,origin = "1970-1-1")))
  return(df)
}

shift <- function(x, lag) {
  n <- length(x)
  xnew <- rep(NA, n)
  if (lag < 0) {
    xnew[1:(n-abs(lag))] <- x[(abs(lag)+1):n]
  } else if (lag > 0) {
    xnew[(lag+1):n] <- x[1:(n-lag)]
  } else {
    xnew <- x
  }
  return(xnew)
}

calculate_titer_change <- function(id, sero_data = as.data.frame(serology), strain = "sH3"){
  cat("id is ", id, "\n")
  data <- sero_data %>% filter(name == strain & memberID == id) %>% 
    select(-c(hhid,member, visit_id))%>%
    group_by(date,name,memberID) %>% 
    summarize(titer = log_titer_trans(max(value))) %>% 
    as.data.frame %>% 
    mutate( date = as.Date(date, origin = '1970-1-1')) %>% 
    arrange(date) 
  if(nrow(data) > 1){
    data$diff <- c(NA,diff(data$titer))
    data$titer_init <- shift(data$titer,1)
    data$delta_t <- c(NA, diff(data$date))
    data$date_midpoint <- data$date - data$delta_t/2
    return(data)
  }
}
