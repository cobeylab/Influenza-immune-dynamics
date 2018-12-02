########################################################################
## Sylvia Ranjeva, August 2018
## Global search in likelihood space 
########################################################################
require(ggplot2)
theme_set(theme_bw(base_size = 18))
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_functions_hh.R")
require(pomp)
require(panelPomp)
select <- dplyr::select
summarize <- dplyr::summarize
rename <- dplyr::rename
n_strains = 1

# IF running on a high-performance computing cluster --------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
chainId = as.numeric(args[1])

## Specify filenames for model input and output files ------------------------------------------------------------------------------------------------------------
output_filename <- paste0("") # Name of the file (".csv") that will store the parameter output and the likelihood
chain_filename <- paste0("chain_",chainId,".rda") # Name of the file (".rda") that will store the MIF chain 
pomp_filename <- "panel_object_H1_all.rda" # Name of the file (".rda") that stores the POMP object

# Specify MIF parameters for this chain ------------------------------------------------------------------------------------------------------------
n_mif = 1
n_mif_updated <- n_mif
n_particles = 2
cooling_rate = .9
n_reps_pfilter = 2
n_particles_pfilter = 2
evaluate_Lhood = TRUE

# Read in and format contact matrix 
source("contact_matrix.R")

## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  include_imprinting = 0,
  imprinting_group = 2,
  measurement_error = 2,
  log_transform_obs = 1,
  include_general_immunity = 1,
  variable_boosting = 1,
  include_k = 1,
  n_strains = 1,
  log_beta_scaled = NA,
  n_age_categories = nrow(contacts),
  age_thres = 15,
  alpha_kids = NA,
  alpha_adults = NA,
  phi_kids = 2.1,
  phi_adults = 2.1,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28), # duration of time between initial inection and peak ab titer,6.8 weeks Zhao et al 2016
  log_k_kids = NA,
  log_k_adults = NA,
  log_w = log(log(2)/204),
  log_r = log(log(2)/(28/10)),
  log_mean_boost_kids = NA,
  log_mean_boost_adults = NA,
  log_sd_boost_kids = NA,
  log_sd_boost_adults = NA,
  logit_d2_kids = NA,
  logit_d2_adults = NA,
  logit_weight_adults = NA,
  logit_weight_kids = NA,
  logit_d_general = NA,
  log_w_general_kids = NA,
  log_w_general_adults = NA,
  log_imprinting_effect_group_1 = NA,
  log_imprinting_effect_group_2 = NA,
  log_omega_hh = NA,
  log_sig = log(1.1),
  log_sig_2 = log(0.74)
)
shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total

rw_sd_vec <- rw.sd(
  log_beta_scaled = 0.01,
  alpha_kids = 0.01,
  alpha_adults = 0.01,
  log_imprinting_effect_group_1 = 0,
  log_imprinting_effect_group_2 = 0,
  logit_d2_adults = 0,
  logit_d2_kids = 0.01,
  logit_d_general = 0,
  log_w_general_kids = 0,
  log_w_general_adults = 0.01,
  logit_weight_adults = 0.01,
  logit_weight_kids = 0,
  log_omega_hh = 0.01
)

inferred_params <- c(
  log_beta_scaled = log(runif(1,0,1)),
  alpha_kids = runif(1,0,5),
  alpha_adults = runif(1,0,5),
  log_k_kids = -.36,
  log_k_adults = -1.1,
  log_imprinting_effect_group_1  = 0,
  log_imprinting_effect_group_2 = 0,#log(runif(1,0.1,5)),
  log_mean_boost_kids = 4.45,
  log_sd_boost_kids = -0.02,
  log_mean_boost_adults = 2.63,
  log_sd_boost_adults = 0.6,
  logit_d2_kids = logit(runif(1,0.1,0.9)),
  logit_d2_adults = -2.5,
  logit_weight_kids = 15,
  logit_weight_adults = logit(runif(1,0.1,0.9)),
  logit_d_general = -50,#logit(runif(1,.1,.9)),
  log_w_general_kids = 100,
  log_w_general_adults = log(runif(1,0,1/180)),
  log_omega_hh = log(runif(1,0,1))
)

source("global_search_methods.R")
 

