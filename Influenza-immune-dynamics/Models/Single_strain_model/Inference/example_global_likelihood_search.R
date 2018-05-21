########################################################################
## Sylvia Ranjeva, May 2018
## Global search in likelihood space 
########################################################################
require(ggplot2)
theme_set(theme_bw(base_size = 18))
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

## Specify filenames for model input and output files ------------------------------------------------------------------------------------------------------------
output_filename <- paste0("") # Name of the file (".csv") that will store the parameter output and the likelihood
chain_filename <- paste0("chain_",chainId,".rda") # Name of the file (".rda") that will store the MIF chain 
pomp_filename <- "panel_object_H1_all.rda" # Name of the file (".rda") that stores the POMP object

# Specify MIF parameters for this chain ------------------------------------------------------------------------------------------------------------
n_mif = 10
n_mif_updated <- n_mif
n_particles = 10e3
cooling_rate = .9
n_reps_pfilter = 3
n_particles_pfilter = 10e3
evaluate_Lhood = TRUE

contacts<- data.frame(age_participant = c(0,10,20,40,65),
                      total = c(1,.953,.746,.751,.562)*7.65)

## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  measurement_error = 2,
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
  log_T_peak = log(28),
  log_w = log(log(2)/204),
  log_r = log(log(2)/(28/10)),
  log_k_1 = NA,  # Antibody ceiling effect for children
  log_k_2 = NA,  # Antibody ceiling effect for adults
  log_imprinting_effect_group_1 = NA, # Fit this parameter if testing imprinting models
  log_imprinting_effect_group_2 = NA, # Fit this parameter if testing imprinting models 
  logit_d2_1 = NA, 
  logit_d2_2 = NA,
  log_mean_boost_1 = NA, # Log mean of short-term titer boost, children
  log_sd_boost_1 = NA, # Log sd of short-term titer boost, children
  log_mean_boost_2 = NA, # Log mean of short-term titer boost, adults
  log_sd_boost_2 = NA, # Log sd of short-term titer boost, adults
  logit_d_general = NA,
  log_w_general_1 = NA,
  log_w_general_2 = NA,
  log_sig = log(1.2),
  log_sig_2 = log(0.74)
)
shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total

rw_sd_vec <- rw.sd(
  log_beta_scaled = 0.01,
  alpha_1 = 0.01,
  alpha_2 = 0.01,
  log_imprinting_effect_group_1 = 0,
  log_imprinting_effect_group_2 = 0,
  logit_d2_1 = 0.01,
  logit_d2_2 = 0.01,
  logit_d_general = 0,
  log_w_general_1 = 0.01,
  log_w_general_2 = 0.01
)

inferred_params <- c(
  log_beta_scaled = log(runif(1,0,1)),
  alpha_1 = runif(1,0,5),
  alpha_2 = runif(1,0,5),
  log_imprinting_effect_group_1  = 0,
  log_imprinting_effect_group_2 = 0,
  logit_d2_1 = logit(runif(1,0.1,0.9)),
  logit_d2_2 = logit(runif(1,0.1,0.9)),
  logit_d_general = -100,
  log_w_general_1 = log(runif(1,0,1/180)),
  log_w_general_2 = log(runif(1,0,1/180))
)

source("global_search_methods.R")
 

