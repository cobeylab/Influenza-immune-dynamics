require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_functions_hh.R")
source("../../../Imprinting/imprinting_functions.R")
require(panelPomp)
select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains 

## Set up input and output files  ## ----------------------------------------------------------------------------------------------------------------------------------
dbFilename <- "../../../Data/Data.sqlite"

## Specify high-level fixed params -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_strains = 1
timestep = 2.5
n_years_initial = 7 # Choose the number of years prior to start of observation period from which to draw last time of infection 
age_threshold = 15 # Division between children and adults
test_subtype = "pH1N1"
calculate_imprinting = T
imprinting_data_filename = NA # Calculate imprinting probabilities and store file (then can set calculate_imprinting = F)
pomp_filename <- paste0("panel_object_", test_subtype,".rda")
community_intensity_table_name <- paste0("community_intensity_", test_subtype)
host_data_filename = paste0("pomp_data_",test_subtype, ".rda")

## Load in contact matrix ---------------------------------------------------------------------------
source("contact_matrix.R")

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Get community intensity data 
db <- dbConnect(SQLite(), dbFilename)
intensity_data <- dbReadTable(db, community_intensity_table_name) %>% 
  mutate(full_date = as.Date(full_date, origin = "1970-1-1"))
dbDisconnect(db)

## Load in the file containing the host data list
load(host_data_filename)

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
  log_beta_scaled = -2.5,
  n_age_categories = nrow(contacts),
  age_thres = age_threshold,
  alpha_kids = log(1000), # wrong value
  alpha_adults = log(1000),
  phi_kids = 2.102,
  phi_adults = 2.102,
  log_mean_gam = log(20),
  log_var_gam = log(1),
  log_T_peak = log(28), # duration of time between initial inection and peak ab titer,6.8 weeks Zhao et al 2016
  log_k_kids = log(.9),
  log_k_adults = log(.9),
  log_w = log(1/180),
  log_r = -1,
  log_mean_boost_kids = log(2),
  log_mean_boost_adults = log(2),
  log_sd_boost_kids = -1000,
  log_sd_boost_adults = -1000,
  logit_d2_kids = -1000,
  logit_d2_adults = -1000,
  logit_d_general = -1000,
  logit_weight_adults = -1000,
  logit_weight_kids = -1000,
  log_w_general_kids = -6,
  log_w_general_adults = -6.5,
  log_imprinting_effect_group_1 = log(1),
  log_imprinting_effect_group_2 = log(1),
  log_omega_hh = .5,
  log_sig = log(1.2),
  log_sig_2 = log(0.74)
)
shared_params[ paste0("age_contact_group_",c(1:nrow(contacts)))] = contacts$age_participant
shared_params[paste0("beta_community_",c(1:nrow(contacts)))] = contacts$total

# Make sure that individuals have more than one visit for longitudinal data 
for( i in c(1:length(host_data_list))){
  if(length(unique(host_data_list[[i]]$dates)) <2){
    host_data_list[[i]] <- NA
  }
  if(length(unique(host_data_list[[i]]$starting_dates$date)) > 1){
    host_data_list[[i]] <- NA
  }
}
ind <- which(is.na(host_data_list))
if(length(ind) > 0 ){
  host_data_list <- host_data_list[-which(is.na(host_data_list))]
}

## Generate pomp object ## --------------------------------------------------------------------------------------
hhs <- c(1:length(host_data_list))

host_data_list <- host_data_list[hhs]
pomp_object_list_complete <-lapply(hhs, 
                           make_pomp_panel_hh, 
                           hh_data_list = host_data_list, 
                           test_params = shared_params, 
                           timestep = timestep, 
                           log_transform_titer = F,
                           n_years_prior = n_years_initial,
                           L_tol = .00001,
                           MAX_MEMBERS = 5,
                           subtype = "pH1N1")
spec_params_df <-data.frame()
n_hh <- length(hhs)
for( i in c(1:n_hh)){
  df <- data.frame(ind = i,
                   param = names(pomp_object_list_complete[[i]]$spec_params),
                   value = as.numeric(pomp_object_list_complete[[i]]$spec_params))
 spec_params_df <- rbind(spec_params_df,df)
  
}

spec_params_input <- reshape(spec_params_df, idvar = "param", timevar= "ind", direction = "wide") 
spec_param_names <- as.character(spec_params_input$param)


if(!file.exists(pomp_filename)){
  pomp_object_list <- lapply(pomp_object_list_complete, function(l) l[[1]])
  spec_params_input = spec_params_input
  spec_param_names <- as.character(spec_params_input$param)
  
  names(pomp_object_list) <- paste0("ind_",c(1:n_hh))
  hh_1 <- pomp_object_list[[1]]
  shared_params <- coef(hh_1)[!(names(coef(hh_1)) %in% spec_param_names )]
  specific_params <- coef(hh_1)[names(coef(hh_1)) %in% spec_param_names]
  spec_params_input <- as.matrix(spec_params_input%>% select(-param))
  rownames(spec_params_input) <- spec_param_names
  colnames(spec_params_input)<- names(pomp_object_list)
  
  ## make panel pomp object from individual objects
  panelPomp(
    object = pomp_object_list,
    shared = shared_params,
    specific = spec_params_input
  ) -> panelObject 

  save(panelObject, host_data_list, file = pomp_filename) 
}
