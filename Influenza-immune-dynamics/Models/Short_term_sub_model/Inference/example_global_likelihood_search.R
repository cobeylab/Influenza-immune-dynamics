########################################################################
## Sylvia Ranjeva May 2018
## Global search in likelihood space 
########################################################################
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_simulation_functions_boosting.R")
require(pomp)
require(panelPomp)
select <- dplyr::select
summarize <- dplyr::summarize
rename <- dplyr::rename
n_strains = 1

## Assign input and output filenames # -------------------------------------------------------------------------------
output_filename <- paste0("global_search_H1_kids.csv") #".csv" file to store search results (params and likelihood)
chain_filename <- paste0("chain_",chainId,"_kids.rda") #".rda" file that stores each MIF chain from a different starting location
pomp_filename <- "panel_object_boosting_H1_kids.rda" #".rda" file containing the pomp object

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

## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  measurement_error = 2,
  log_transform_obs = 1,
  variable_boosting = 1,
  include_k = 1,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28), 
  log_k = NA,
  log_w = log(log(2)/204),
  log_r = log(log(2)/(28/10)),
  log_mean_boost = NA,
  log_sd_boost = NA,
  log_sig = NA,
  log_sig_2 = NA
)


rw_sd_vec <- rw.sd(
  log_mean_boost = 0.01,
  log_sd_boost = 0.01,
  log_k = 0.01
)

inferred_params <- c(
 log_mean_boost = log(runif(1,0,100)),
 log_sd_boost = log(runif(1,0,100)),
 log_k = log(runif(1,0,10)),
 log_sig = log(.6),
 log_sig_2 = log(.74)
)

source("global_search_methods.R")
 

