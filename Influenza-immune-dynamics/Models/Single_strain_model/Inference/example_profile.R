########################################################################
## Sylvia Ranjeva, May 2018
## Example likelihood profile over one parameter  
########################################################################
require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_functions.R")
require(pomp)
require(panelPomp)
select <- dplyr::select
summarize <- dplyr::summarize
rename <- dplyr::rename
n_strains = 1

# IF running on a high-performance computing cluster --------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
chainId = as.numeric(args[1])
#------------------------------------------------------------------------------------------------------------
# Specify output filenames
output_filename <- "" # Name of file (".csv" file) to store end parameters and likelihood
chain_filename <- paste0("chain_",chainId,".rda") # Name of file (".rda" file) to store MIF chain
pomp_filename <-  "" # Name of file (".rda" file) containing the POMP object

# Specify MIF parameters for this chain
n_mif = 100
n_mif_updated <- n_mif
n_particles = 10e3
cooling_rate = .9
n_reps_pfilter = 10 # Number of replicate passes through the particle filter after MIF search
n_particles_pfilter = 20e3
evaluate_Lhood = TRUE
n_profile_points = 100
n_profile_reps = 3 # Number of replicate points for each value of the focal parameter
alpha_min = -4 # Range over which to profile 
alpha_max = -2 # Range over which to profile 
alpha_vals = rep(seq(alpha_min, alpha_max, length.out = n_profile_points), n_profile_reps)

# Contact matrix
contacts<- data.frame(age_participant = c(0,10,20,40,65),
                      total = c(1,.953,.746,.751,.562)*7.65)
## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  include_imprinting = 0,
  imprinting_group = 2,
  log_transform_obs = 1,
  include_general_immunity = 1,
  include_k = 1,
  variable_boosting = 1,
  n_strains = 1,
  log_beta_scaled = NA,
  n_age_categories = nrow(contacts),
  age_thres = 15,
  alpha_1 = NA, 
  alpha_2 = NA,
  phi_1 = 2.1,
  phi_2 = 2.1,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28), # duration of time between initial inection and peak ab titer, Zhao 2016
  log_w = log(log(2)/204),
  log_r = log(log(2)/(28/10)),
  log_k_1 = NA,
  log_k_2 = NA,
  log_imprinting_effect_group_1 = NA,
  log_imprinting_effect_group_2 = NA,
  logit_d2_1 = NA,
  logit_d2_2 = NA,
  log_mean_boost_1 = NA,
  log_sd_boost_1 = NA,
  log_mean_boost_2 = NA,
  log_sd_boost_2 = NA,
  logit_d_general = NA,
  log_w_general_1 = NA,
  log_w_general_2 = NA,
  log_sig = log(1.29),
  log_sig_2 = log(0.74)
)
shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total

rw_sd_vec <- rw.sd(
  log_beta_scaled = 0, # Example profile over beta (transmission rate) shown here. Set rw.sd of focal param to zero.
  alpha_1 = 0.01,
  alpha_2 = 0.01,
  log_k_1 = 0,
  log_k_2 = 0,
  log_imprinting_effect_group_1 = 0,
  log_imprinting_effect_group_2 = 0,
  log_mean_boost_1 = 0,
  log_sd_boost_1 = 0,
  log_mean_boost_2 = 0,
  log_sd_boost_2 = 0,
  logit_d2_1 = 0.01,
  logit_d2_2 = 0.01,
  logit_d_general = 0,
  log_w_general_1 = 0.01,
  log_w_general_2 = 0.01
)

inferred_params <- c(
  log_beta_scaled = alpha_vals[chainId], # Example profile over beta (transmission rate) shown here. Fix focal param.
  alpha_1 = runif(1,0,5),
  alpha_2 = runif(1,0,5),
  log_k_1 = -0.5,
  log_k_2 = -1.2,
  log_imprinting_effect_group_1  = 0,
  log_imprinting_effect_group_2 = 0,#log(runif(1,0.1,5)),
  log_mean_boost_1 = 3.5,
  log_sd_boost_1 = -0.1,
  log_mean_boost_2 = 3.6,
  log_sd_boost_2 = 0.6,
  logit_d2_1 = logit(runif(1,0.1,0.9)),
  logit_d2_2 =  logit(runif(1,0.1,0.9)),
  logit_d_general = -50,#logit(runif(1,.1,.9)),
  log_w_general_1 = log(runif(1,0,1/180)),
  log_w_general_2 = log(runif(1,0,1/180))
)

source("global_search_methods.R")
 

