
## ----------------------------------------------------------------------------------
## Sylvia Ranjeva, August 2018
## Simulate latent state data from the best-fit parameter set.
## ----------------------------------------------------------------------------------

require(ggplot2)
theme_set(theme_bw(base_size = 18))
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../../Utility_scripts/model_functions_hh.R")
require(pomp)
require(panelPomp)
select = dplyr::select
summarize = dplyr::summarize
rename = dplyr::rename
group_by = dplyr::group_by

# IF running on a high-performance computing cluster --------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
chainId = as.numeric(args[1])

## Input params/specs, filenames, etc. ## ---------------------------------------------------------------------------------------------------------------------
n_strains = 1
pomp_object_filename <- "" # Name of file (".rda" file) containing pomp object
n_sim = 1000 # Desired number of simulations
this_subtype = "pH1N1"
simulation_output_filename <- paste0("simulated_data_",n_sim, "_sim_",this_subtype,".rda")
simulation_params_filename <- "" # Name of file (".rda" file) containing simulation parameters, which should be a named vector called "pars"

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
load(pomp_object_filename)
load(simulation_params_filename)

titer_df_all <- data.frame()
inf_times_df_all <- data.frame()
for( j in c(1:n_sim)){
  cat("j is : ", j, "\n")
  for( i in c(1:length(panelObject))){
    hh1 <- panelObject@unit.objects[[i]]
    dat <- host_data_list[[i]]
    hh_params <- hh1@params[!(names(hh1@params) %in% names(panelObject@pParams$shared))]
    n_members = as.numeric(hh_params["n_members"])
    members_hh = dat$demog_init$member
    s1 <- simulate(hh1, params = c(pars,hh_params))
    df_s1 <- as.data.frame(t(states(s1)))
    for( k in c(1:n_members)){
      df_ind <- df_s1 %>% select(contains(paste0("m",k)))
      demog_ind <- dat$demog_init %>% filter(member == members_hh[k])
      titer_df <- data.frame(hh = i,
                             member = k,
                             sim = j,
                             time = s1@times,
                             date = as.Date(as.Date(as.character(dat$dates[1]),origin = "1970-1-1") + as.numeric(hh1@times), origin = '1970-1-1'),
                             age = demog_ind$age,
                             h_latent = sapply(df_ind$h_m,log_titer_trans), 
                             h_baseline = sapply(df_ind$h_baseline_individual,log_titer_trans),
                             q_latent = df_ind$q_m
                             )
      titer_df_all <- rbind(titer_df_all, titer_df)
      
      date_time_df <- data.frame(time = as.numeric(hh1@times),
                                 date = as.Date(as.Date(as.character(dat$dates[1]),origin = "1970-1-1") + as.numeric(hh1@times), origin = '1970-1-1')
      )
      
      inf_times_sim = unique(df_ind$t_last_infection_m) 
      n_inf <- sum(inf_times_sim > -1)
      
      if(n_inf > 0){
        inf_times_df <- data.frame( hh = i,
                                    member = k,
                                    age = demog_ind$age,
                                    inf_time = inf_times_sim[inf_times_sim > -1] 
        )
        inf_times_df$date = as.Date(as.numeric(as.Date(date_time_df$date[1], origin = '1970-1-1')) + inf_times_df$inf_time, origin = '1970-1-1') 
        inf_times_df$sim = j
        inf_times_df_all <- rbind(inf_times_df_all, inf_times_df)
      }
   }
  }
}


titer_df_all <-  titer_df_all %>%  mutate(ind = group_indices_(., .dots=c("hh", "member")))
inf_times_df_all$ind = NA
for(i in c(1:nrow(inf_times_df_all))){
  inf_times_df_all[i,]$ind = unique(titer_df_all[titer_df_all$hh == inf_times_df_all[i,]$hh & titer_df_all$member == inf_times_df_all[i,]$member,]$ind)
} 

save(inf_times_df_all, titer_df_all, pars, file = simulation_output_filename)
