require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_functions.R")
require(panelPomp)
select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains 

## Set up input and output files  ## ----------------------------------------------------------------------------------------------------------------------------------
test_data_filename <- "" # File (".rda" file) containing relevant data 
pomp_filename <- " " # Name of file (".rda" file) storing pomp object
population_data_filename <- ".././Data/L_data_H1.rda"

## Specify high-level fixed params -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
contacts<- data.frame(age_participant = c(0,10,20,40,65),
                      total = c(1,.926,.847,.812,.608)*12.62)
n_strains = 1
timestep = 2.5
n_years_initial = 7 # Choose the number of years prior to start of observation period from which to draw last time of infection 
age_threshold = 15 # Division between children and adults
subtype = "pH1N1"
age_group = "children" # Choose between "children" and "adults"
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

load(test_data_filename)
#Load in population level data 
load(population_data_filename)

host_ages <- data.frame() # Specify host ages so that we can filter children vs adults
for( i in c(1:length(host_data_list))){
    df_ind <- data.frame( ind = i, age = host_data_list[[i]]$demog_init["age"])
    host_ages <- rbind(imp_probs,df_ind)
}

# Choose individuals for this analysis based on age group 
if(age_group == "children" ){
  individuals <- which(host_ages$age <= age_threshold)
}

host_data_list <- host_data_list[individuals]

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
  log_beta_scaled = log(5),
  n_age_categories = nrow(contacts),
  age_thres = 19,
  alpha_1 = log(40), # wrong value
  alpha_2 = log(40),
  phi_1 = 2.102,
  phi_2 = 2.192,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28),
  log_k_1 = log(.9),
  log_k_2 = log(.9),
  log_w = log(1/180),
  log_r = -1,
  log_mean_boost_1 = log(8),
  log_mean_boost_2 = log(8),
  log_sd_boost_1 = log(4),
  log_sd_boost_2 = log(4),
  logit_d2_1 = logit(.5),
  logit_d2_2 = logit(.5),
  logit_d_general = logit(.5),
  log_w_general_1 = log(1/50),
  log_w_general_2 = log(1/50),
  log_imprinting_effect_group_1 = log(1),
  log_imprinting_effect_group_2 = log(1),
  log_D = log(1e10),
  log_sig = log(.5),
  log_sig_2 = log(1)
)
shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total

# Make sure that individuals have more than one visit for longitudinal data 
for( i in c(1:length(host_data_list))){
  if(length(unique(host_data_list[[i]]$dates)) <2){
    host_data_list[[i]] <- NA
  }
}
ind <- which(is.na(host_data_list))
if(length(ind) > 0 ){
  host_data_list <- host_data_list[-which(is.na(host_data_list))]
}


pomp_object_list <-lapply(c(1:length(host_data_list)), 
                           make_pomp_panel, 
                           ind_data_list = host_data_list, 
                           test_params = shared_params, 
                           timestep = timestep, 
                           log_transform_titer = F,
                           n_years_prior = n_years_initial,
                           L_tol = .00001,
                           this_subtype = subtype)
spec_params_df <-data.frame()
for( i in c(1:n_ind)){
  df <- data.frame(ind = i,
                   param = names(pomp_object_list[[i]]$spec_params),
                   value = as.numeric(pomp_object_list[[i]]$spec_params))
 spec_params_df <- rbind(spec_params_df,df)
  
}

spec_params_input <- reshape(spec_params_df, idvar = "param", timevar= "ind", direction = "wide") 
spec_param_names <- as.character(spec_params_input$param)

n_ind <- length(inds)
if(!file.exists(pomp_filename)){
  pomp_object_list <- lapply(pomp_object_list, function(l) l[[1]])
  names(pomp_object_list) <- paste0("ind_",c(1:n_ind))
  hh_1 <- pomp_object_list[[1]]
  shared_params <- coef(hh_1)[!(names(coef(hh_1)) %in% spec_param_names )]
  specific_params <- coef(hh_1)[names(coef(hh_1)) %in% spec_param_names]
  spec_params_input <- as.matrix(spec_params_input%>% select(-param))
  rownames(spec_params_input) <- spec_param_names
  colnames(spec_params_input)<- names(pomp_object_list)
  
  ## Make panel pomp object from individual objects
  panelPomp(
    object = pomp_object_list,
    shared = shared_params,
    specific = spec_params_input
  ) -> panelObject 

  save(panelObject, host_data_list, file = pomp_filename) 
}
