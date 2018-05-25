
## ----------------------------------------------------------------------------------
## Sylvia Ranjeva, May 2018
## Simulate data (separately for kids and adults) from the best-fit parameter set.
## ----------------------------------------------------------------------------------

require(ggplot2)
theme_set(theme_bw(base_size = 18))
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../../Utility_scripts/model_functions.R")
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
pomp_object_filename <- "../../Inference/FILEAME" # Name of file (".rda" file) containing pomp object. This file is stored in the Inference folder
n_sim = 1000 # Desired number of simulations
age_group = "adults" # Select "children" or "adults"
this_subtype = "pH1N1"
simulation_output_filename <- paste0("simulated_data_",n_sim, "_validation_",this_subtype,"_", age_group,".rda")
simulation_params_filename <- "" # Name of file (".rda" file) containing simulation parameters, which should be a named vector called "pars"
## --------------------------------------------------------------------------------------------------------------------------------------------------------------
load(simulation_params_filename)
load(pomp_object_filename)

titer_df_all <- data.frame()
inf_times_df_all <- data.frame()
for( j in c(1:n_sim)){
  cat("j is : ", j, "\n")
  for( i in c(1:length(panelObject))){
    ind1 <- panelObject@unit.objects[[i]]
    dat <- host_data_list[[i]]
    ind_params <- unlist(ind1@params[c("init_age","h_t0_s1","p_imprinted_h3", "p_imprinted_group_1","h_baseline_individual_observed")])
    s1 <- simulate(ind1, params = c(pars,ind_params ))
    titer_df <- data.frame(ind = i,
                           sim = j,
                           time = s1@times,
                           date = as.Date(as.Date(as.character(dat$dates[1]),origin = "1970-1-1") + as.numeric(ind1@times), origin = '1970-1-1'),
                           age = dat$demog_init[1],
                           h_latent = sapply(as.numeric(states(s1)["h_s1",]),log_titer_trans), 
                           h_baseline = sapply(as.numeric(states(s1)["h_baseline_individual",]),log_titer_trans),
                           q_latent = as.numeric(states(s1)["q",]
                           )
    )
    titer_df_all <- rbind(titer_df_all, titer_df)
    
    date_time_df <- data.frame(time = as.numeric(ind1@times),
                               date = as.Date(as.Date(as.character(dat$dates[1]),origin = "1970-1-1") + as.numeric(ind1@times), origin = '1970-1-1')
    )
    
    inf_times_sim = unique(as.numeric(states(s1)["t_last_infection_s1",])) 
    n_inf <- sum(inf_times_sim > -1)
    
    if(n_inf > 0){
      inf_times_df <- data.frame( ind = i,
                                  age = dat$demog_init[1],
                                  inf_time = unique(as.numeric(states(s1)["t_last_infection_s1",]))) %>% filter(inf_time > -1)
      
      inf_times_df$date = as.Date(as.numeric(as.Date(date_time_df$date[1], origin = '1970-1-1')) + inf_times_df$inf_time, origin = '1970-1-1') 
      inf_times_df$sim = j
      inf_times_df_all <- rbind(inf_times_df_all, inf_times_df)
    }
  }
}

save(inf_times_df_all, titer_df_all, pars, file = simulation_output_filename)
