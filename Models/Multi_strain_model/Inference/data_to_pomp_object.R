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

dbFilename <- "../../../Data/Data.sqlite"
pomp_filename <- "panel_object_H1_H3.rda" # Output file (".rda") storing pomp object
test_data_filename <- "pomp_data_H1_H3.rda" # Input file (".rda") storing host data list
timestep = 2.5 # Timestep for the simulation, days
n_strains = 2

## Load in host data list
load(test_data_filename)

## Get community intensity data 
db <- dbConnect(SQLite(), dbFilename)
intensity_data_H3 <- dbReadTable(db, "community_intensity_H3N2") %>% 
  mutate(full_date = as.Date(full_date, origin = "1970-1-1"))
intensity_data_H1 <- dbReadTable(db, "community_intensity_pH1N1") %>% 
  mutate(full_date = as.Date(full_date, origin = "1970-1-1"))
dbDisconnect(db)
intensity_data <- rbind(intensity_data_H1, intensity_data_H3)

for( i in c(1:length(host_data_list))){
  if(length(unique(host_data_list[[i]]$dates)) <2){
    host_data_list[[i]] <- NA
  }
}
ind <- which(is.na(host_data_list))
if(length(ind) > 0 ){
  host_data_list <- host_data_list[-which(is.na(host_data_list))]
}


## Generate test params 
contacts<- data.frame(age_participant = c(0,10,20,40,65),
                      total = c(1,.953,.746,.751,.562)*7.65)
## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  include_imprinting = 0,
  include_heterosubtypic = 1,
  log_transform_obs = 1,
  include_general_immunity = 0,
  include_k = 1,
  variable_boosting = 1,
  n_strains = 1,
  n_age_categories = nrow(contacts),
  age_thres = 25,
  phi_1 = 2.1,
  phi_2 = 2.1,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28), 
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

for( i in c(1:length(host_data_list))){
  if(length(unique(host_data_list[[i]]$dates)) <2){
    host_data_list[[i]] <- NA
  }
}
ind <- which(is.na(host_data_list))
if(length(ind) > 0 ){
  host_data_list <- host_data_list[-which(is.na(host_data_list))]
}

n_ind <- length(host_data_list)
inds <- c(1:n_ind)

pomp_object_list <- lapply(inds, 
                           make_pomp_panel_2strain, 
                           ind_data_list = host_data_list, 
                           test_params = shared_params, 
                           timestep = timestep, 
                           log_transform_titer = F,
                           L_data = intensity_data,
                           n_years_prior = 7,
                           L_tol = .00001,
                           subtypes = c("pH1N1","H3N2")
                           )
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
  
  ## make panel pomp object from individual objects
  panelPomp(
    object = pomp_object_list,
    shared = shared_params,
    specific = spec_params_input
  ) -> panelObject 

  save(panelObject, host_data_list,  file = pomp_filename) 
}
