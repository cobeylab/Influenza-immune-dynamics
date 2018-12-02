
log_titer_trans <- function(x){
  return(log2(x/10) + 2)
}

get_age_at_recruitment <- function(memberID, demography_table = demography){
  return(unique(demography_table[demography_table$memberID == memberID,]$age_at_recruitment))
}

get_birth_date <- function(memberID, demography_table = demography){
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

get_initial_age_for_ind <- function(member_num, hh_num, dem_data = demography, DAYSPERYEAR = 365){
  data <- dem_data %>% filter(hhid == hh_num & member == member_num)
  dem_df <- data.frame(member = data$member,
                       age = as.numeric(as.Date(data$visitdate) - as.Date(data$birthdate))/DAYSPERYEAR)
  return(dem_df)
}

get_serology_data_for_ind_and_strain <- function(strain_id, id , strain_vec = strains, serology_data = serology){
  strain_name <- strain_vec[strain_id]
  data <- serology %>% filter(memberID == id & subtype == strain_name) %>% as.data.frame %>% select(-c(hhid,member,visit_id,subtype)) %>% 
    filter(!is.na(date))
  n_strains <- length(strain_vec)
  data <- data %>% select(-memberID) %>% arrange(date)
  names(data) <- c("date",paste0("h_obs_s",strain_id))
  data$obs_time <- (as.numeric(as.Date(data$date,origin = '1970-1-1')) - as.numeric(as.Date(data$date[1],origin = '1970-1-1')))
  data <- data %>% arrange(obs_time)
  data <- melt(data, id.vars = c("obs_time","date"))
  return(data)
}

get_hh_serology_data_for_ind_and_strain <- function(member_num, id, strain_id = 1, strain_vec = strains, serology_data = serology, member_map = member_df ){
  strain_name <- strain_vec[strain_id]
  data <- serology_data %>% filter(hhid == id & member == member_num & subtype == strain_name) %>% as.data.frame %>% select(-c(hhid,visit_id,subtype)) %>% 
    filter(!is.na(date))
  n_strains <- length(strain_vec)
  data <- data %>% select(-memberID) %>% arrange(date)
  mem_num <- member_map[which(member_map$mem_orig == member_num),]$mem_new
  names(data) <- c("member","date", paste0("h_obs_m",mem_num))
  data$obs_time <- (as.numeric(as.Date(data$date,origin = '1970-1-1')) - as.numeric(as.Date(data$date[1],origin = '1970-1-1')))
  data <- data %>% arrange(obs_time)
  data <- melt(data, id.vars = c("member","obs_time","date"))
}

## Make pomp data

make_pomp_data_for_hh <- function(hh_num, 
                                   strain_vec = strains, 
                                   serology_data = serology, 
                                   imprinting_data,
                                   simulate_times = F, 
                                   demographic_data = demography,
                                   n_vis_sim = NA){
  cat("hh is ", hh_num, "\n")
  n_strains <- length(strain_vec)
  serology_subset <- serology %>% filter(hhid == hh_num)
  members <- sort(unique(serology_subset$member))
  member_df <- data.frame(mem_orig = members, mem_new = c(1:length(members)))
  hh_data_titer_list <- lapply(members, get_hh_serology_data_for_ind_and_strain, id = hh_num, member_map = member_df)

  hh_data_titer <- do.call("rbind",hh_data_titer_list)
  
  hh <- hh_data_titer %>% 
    group_by(member) %>% 
    mutate(date_diff = c(as.numeric(diff(date)),NA)) %>% 
    filter(date_diff > 20 | is.na(date_diff)) %>% 
    group_by(variable) %>% 
    mutate(obs_time = obs_time - min(obs_time))
  
  starting_dates <- hh %>% filter(obs_time == 0) %>% as.data.frame() %>% select(c(member,date))
  
  data_all <- hh %>% select(-date_diff) %>% as.data.frame()
  dates <- unique(data_all$date)
  data <- reshape(data_all %>% select(-c(date,member)), timevar = "variable", idvar = "obs_time", direction = "wide" )
  names(data) <- c("obs_time", as.character(unique(data_all$variable)))
  #data <- data %>% drop_na()
  demog_init <- do.call("rbind",lapply(members, get_initial_age_for_ind, hh_num = hh_num, dem_data = demographic_data))
  imprinting_probs <- imprinting_data %>% filter(memberID %in% serology_subset$memberID) %>% 
    mutate(member = NA)
  for( i in c(1:nrow(imprinting_probs))){
    imprinting_probs[i,]$member = unique(serology_subset[which(serology_subset$memberID == imprinting_probs[i,]$memberID),]$member)
  }
  
  min_titer <- data_all %>% 
    group_by(member) %>% 
    summarize(min = min(value))
  
   if(simulate_times){
    init_titer <- hh_data_titer %>% filter(obs_time == 0)
    date_0 <- dates[1]
    t_0 <- min(data$obs_time)
    t_end <- max(data$obs_time)
    obs_time <-  seq(t_0, t_end, length.out = n_vis_sim)
    data = data.frame(obs_time = obs_time)
    for(i in c(1:length(members))){
      min_val = min_titer[min_titer$member == members[i],]$min
      init_titer_val <- init_titer %>% filter(init_titer$variable == paste0("h_obs_m",i)) %>% select(value) %>% as.numeric
      data  <- cbind(data,c(init_titer_val,rep(min_val, n_vis_sim -1)))
    }
    names(data) <- c("obs_time",paste0("h_obs_m",c(1:length(unique(members)))))
    dates <- as.Date(as.numeric(as.Date(date_0, origin ='1970-1-1')) + obs_time, origin = '1970-1-1')
  }
  return(list(y = data, n_strains = n_strains, n_members = length(unique(members)), demog_init = demog_init, dates = dates, starting_dates = starting_dates, imprinting_probs = imprinting_probs))
} 



make_pomp_data_for_ind <- function(id, 
                                   strain_vec = strains, 
                                   serology_data = serology, 
                                   imprinting_data,
                                   simulate_times = F, 
                                   demographic_data = demography,
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
  data <- data %>% drop_na()
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
