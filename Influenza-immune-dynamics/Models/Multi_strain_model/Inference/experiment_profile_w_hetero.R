###
#####################################################################
## Global search in likelihood space 
########################################################################
require(ggplot2)
theme_set(theme_bw(base_size = 18))
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source(".././Utility_scripts/model_functions.R")
require(pomp)
require(panelPomp)
select <- dplyr::select
summarize <- dplyr::summarize
rename <- dplyr::rename
n_strains = 2

# IF running on a high-performance computing cluster --------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
chainId = as.numeric(args[1])
#------------------------------------------------------------------------------------------------------------

# Specify MIF parameters for this chain
n_mif = 100
n_mif_updated <- n_mif
n_particles = 10e3
cooling_rate = .9
n_reps_pfilter = 3
n_particles_pfilter = 20e3
evaluate_Lhood = TRUE
n_profile_points = 100
n_profile_reps = 1
alpha_min = -2
alpha_max = -.5
alpha_vals = rep(seq(alpha_min, alpha_max, length.out = n_profile_points), n_profile_reps)

# Format contact matrix 
contacts<- data.frame(age_participant = c(0,10,20,40,65),
                      total = c(1,.953,.746,.751,.562)*7.65)

## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  include_imprinting = 0,
  include_heterosubtypic = 1,
  log_transform_obs = 1,
  include_general_immunity = 1,
  include_k = 1,
  variable_boosting = 1,
  n_strains = 2,
  n_age_categories = nrow(contacts),
  age_thres = 15,
  phi_1 = 2.1,
  phi_2 = 2.1,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28), # duration of time between initial inection and peak ab titer, Zhao 2016
  log_w = log(log(2)/204),
  log_r = log(log(2)/(28/10)),
  log_w_hetero = NA,
  logit_d_general = -50,
  log_imprinting_effect_group_1 = log(1),
  log_imprinting_effect_group_2 = log(1),
  log_sig = log(1.29),
  log_sig_2 = log(0.74)
)

shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total

shared_params[paste0("imprinting_group_s",c(1:n_strains))] = c(1,2)
shared_params[paste0("log_beta_scaled_s",c(1:n_strains))] = c(-2.8,-3.1)
shared_params[paste0("log_mean_boost_1_s",c(1:n_strains))] = c(3.5, 5.1) 
shared_params[paste0("log_mean_boost_2_s",c(1:n_strains))] = c(3.6, 4.6) 
shared_params[paste0("log_sd_boost_1_s",c(1:n_strains))] = c(-0.1,0.4)
shared_params[paste0("log_sd_boost_2_s",c(1:n_strains))] = c(0.6, 0.3)
shared_params[paste0("logit_d2_1_s",c(1:n_strains))] = c(-1.33,-3.9)
shared_params[paste0("logit_d2_2_s",c(1:n_strains))] = c(-100,-100) # No evidence of long-term boost in adults
shared_params[paste0("log_w_general_1_s",c(1:n_strains))] = 100 #No nonspecific immunity in kids
shared_params[paste0("log_w_general_2_s",c(1:n_strains))] = c(-7.6, -7.0)
shared_params[paste0("logit_c_s",c(1:n_strains))] = c(NA,NA) 
shared_params[paste0("log_k_1_s",c(1:n_strains))] = c(-.5,-.7)
shared_params[paste0("log_k_2_s",c(1:n_strains))] = c(-1.2,0.0)
shared_params[paste0("alpha_1_s",c(1:n_strains))] = c(4.0,4.75)
shared_params[paste0("alpha_2_s",c(1:n_strains))] = c(5.8,6.3)
shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total


output_filename <- paste0("profile_w_hetero.csv")
chain_filename <- paste0("chain_",chainId,"_prof_w_hetero.rda") 
pomp_filename <-  "panel_object_H1_H3.rda" 

rw_sd_vec <- rw.sd(
  log_w_hetero = 0,
  logit_c_s1 = 0,
  logit_c_s2 =0
)

inferred_params <- c(
  log_w_hetero = alpha_vals[chainId],
  logit_c_s1 = -100, #Don't allow partial immunity
  logit_c_s2 = -100  #Don't allow partial immunity
)
source("global_search_methods.R")
 


