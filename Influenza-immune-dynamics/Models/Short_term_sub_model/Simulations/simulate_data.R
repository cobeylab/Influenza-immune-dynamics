
########################################################################
## Sylvia Ranjeva, May 2018
## Simulate data from the best-fit parameter set for the sub-model
########################################################################
require(ggplot2)
theme_set(theme_bw(base_size = 18))
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_simulation_functions.R")
require(pomp)
require(panelPomp)
select <- dplyr::select

## Input params/specs, filenames, etc. ## ---------------------------------------------------------------------------------------------------------------------
pomp_object_filename <- "" # named pomp object R data file (.rda)
pomp_data_filename <- "" # named pomp data R data file (.rda)
n_strains = 1 # sub-model is currently implemented for single strains
simulation_output_filename <- "" # named R data file (.rda)
parameter_filename <- "" # named R data file (.rda) with simulation outputs
n_sim =  # number of simulations

## Data ## ----------------------------------------------------------------------------------------------------
load(pomp_data_filename)
serology_sub <- ## SLR: code to extract relevatn datea table# 
ids <- demog %>% filter(age_at_recruitment > 15) %>% select(memberID) %>% unlist() %>% as.numeric()
serology_sub <- serology_sub %>% filter(memberID %in% ids)
load(pomp_object_filename)
load(parameter_filename)
## --------------------------------------------------------------------------------------------------------------------------------------------------------------

coef <- coef(panelObject)$shared
for(i in c(1:length(pars_shared))){
  coef[which(names(coef) == names(pars_shared)[i])] <- pars_shared[i]
}


titer_df_all <- data.frame()
inf_times_df_all <- data.frame()
test_ids <- unique(serology_sub$memberID)

for( j in c(1:n_sim)){
  cat("j is : ", j, "\n")
  for( i in c(1:length(panelObject))){
    ind1 <- panelObject@unit.objects[[i]]
    dat <- serology_sub %>% filter(memberID == test_ids[i])
    dem <- demog %>% filter(memberID == test_ids[i])
    ind_params <- unlist(ind1@params[c("init_age","h_t0","t_infection")])
    s1 <- simulate(ind1, params = c(coef,ind_params))
    titer_df <- data.frame(ind = i,
                           sim = j,
                           time = s1@times,
                           t_inf = coef(s1)["t_infection"],
                           date_inf = as.Date(as.character(dat$date[1]),origin = "1970-1-1") + as.numeric(coef(s1)["t_infection"]),
                           date = as.Date(as.Date(as.character(dat$date[1]),origin = "1970-1-1") + as.numeric(ind1@times), origin = '1970-1-1'),
                           age = as.numeric(dem$age_at_recruitment),
                           h_latent = sapply(as.numeric(states(s1)["h",]),log_titer_trans), 
                           h_obs_sim = as.numeric(obs(s1)["h_obs",]),
                           h_obs = sapply(as.numeric(obs(ind1)["h_obs",]),log_titer_trans),
                           boost_obs_sim = c(0, diff(as.numeric(obs(s1)["h_obs",]))),
                           boost_latent = c(0, diff(log_titer_trans(as.numeric(states(s1)["h",])))),
                           boost_obs =  c(0, diff(log_titer_trans(as.numeric(obs(ind1)["h_obs",]))))
                           )
    
    titer_df_all <- rbind(titer_df_all, titer_df)
    date_time_df <- data.frame(time = as.numeric(ind1@times),
                               date = as.Date(as.Date(as.character(dat$date[1]),origin = "1970-1-1") + as.numeric(ind1@times), origin = '1970-1-1')
    )
    
  }
}
save( titer_df_all, file = simulation_output_filename)
