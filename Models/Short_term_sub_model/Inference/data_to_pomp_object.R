require(ggplot2)
theme_set(theme_bw(base_size = 18))
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

## Specify the name of the pomp object file ------------------------------------------------------
test_subtype = "pH1N1" # Choose between "pH1N1" and "H3N2"
data_table_name <- paste0("data_sub_model_",test_subtype)
pomp_filename <- "" # Name the pomp object file (".rda" file)
dbFilename <- "../../../Data/Data.sqlite"
#Read in data  
db <- dbConnect(SQLite(), dbFilename)
serology_sub <- dbReadTable(db, data_table_name) 
dbDisconnect(db)


# Filtering for kids <= 15 y old, but can be altered to generate pomp object for adults too
serology_sub <- serology_sub %>% filter(age_at_recruitment <= 15)

timestep = 2.5 # Simulation timestep, days

## Generate test params 
shared_params <- c(
  log_transform_titers = 1,
  log_transform_obs = 1,
  variable_boosting = 1,
  include_k = 1,
  log_mean_gam = log(5),
  log_var_gam = log(1),
  log_T_peak = log(28), 
  log_k = log(.9), # trivial value just to construct the pomp object
  log_w = log(1/180),
  log_r = log(log(2)/(28*10)),
  log_mean_boost = log(8), # trivial value just to construct the pomp object
  log_sd_boost = log(.1), # trivial value just to construct the pomp object
  log_sig = log(1.29),
  log_sig_2 = log(0.74)
)

test_ids = unique(serology_sub$memberID)

pomp_object_list <- lapply(c(1:length(test_ids)),
                           test_ids = test_ids,
                           make_pomp_panel,  
                           data = serology_sub,
                           test_params = shared_params, 
                           timestep = timestep, 
                           log_transform_titer = F
)
spec_params_df <-data.frame()

for( i in c(1:length(test_ids))){
  df <- data.frame(ind = i,
                   param = names(pomp_object_list[[i]]$spec_params),
                   value = as.numeric(pomp_object_list[[i]]$spec_params))
 spec_params_df <- rbind(spec_params_df,df)
  
}

spec_params_input <- reshape(spec_params_df, idvar = "param", timevar= "ind", direction = "wide") 
spec_param_names <- as.character(spec_params_input$param)
n_ind = length(test_ids)
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

  save(panelObject,  file = pomp_filename) 
}
