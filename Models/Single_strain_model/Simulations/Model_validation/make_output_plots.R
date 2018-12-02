########################################################################
require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(lubridate)
require(dplyr)
source("../../Utility_scripts/model_functions_hh.R")
textSize = 12
source("../../Utility_scripts/plot_themes.R")
require(pomp)
require(panelPomp)
require(cowplot)
select <- dplyr::select
summarize <- dplyr::summarize
mutate <- dplyr::mutate 
group_by <- dplyr::group_by

## Input specs ## ----------------------------------------------------------------------------------------------------------------------------
age_threshold = 15 # Age that divides kids and adults
this_subtype = "pH1N1" # Choose either "pH1N1" or "H3N2" 
n_strains = 1
save_plots <- T # Save plots to pdf?

# Get the community intensity data 
intensity_data_tablename <- paste0("community_intensity_", this_subtype)
dbFilename <- "../../../../Data/Data.sqlite"
db <- dbConnect(SQLite(), dbFilename)
intensity_data <- dbReadTable(db, intensity_data_tablename) %>% 
  mutate(full_date = as.Date(full_date, origin = "1970-1-1"))
dbDisconnect(db)

# Load in the simulated host list and simulated output 
pomp_object_filename <- "../../Inference/FILEAME" # Name of file (".rda" file) containing pomp object. This file is stored in the Inference folder
model_predictions_filename <- "" # Name of file (".rda") containing model simulations 

load(pomp_object_filename)
rm(panelObject)
load(model_predictions_filename)


## Apply measurement error --------------------------------------------------------------------
titer_df_all$h_obs_sim = sapply(titer_df_all$h_latent, err_func, thres = 1, sig = 1.2, sig2 = 0.74)

## Format data ##----------------------------------------------------------------------------------------------------------------------------

for( i in c(1:length(host_data_list))){
  n_vis <- nrow(host_data_list[[i]]$y)
  if(length(unique(host_data_list[[i]]$dates))< 2){
    host_data_list[[i]] <- NA
  }
}
ind <- which(is.na(host_data_list))
if(length(ind) > 0 ){
  host_data_list <- host_data_list[-which(is.na(host_data_list))]
}

age_vec <- array(NA, length(host_data_list))
df_demog = data.frame()
for( i in c(1:length(host_data_list))){
  df_demog <- rbind(df_demog, data.frame(hh= i,
                                         member = c(1:length(host_data_list[[i]]$demog_init$member)),
                                         first_date = host_data_list[[i]]$dates[1],
                                         last_date = host_data_list[[i]]$dates[length(host_data_list[[i]]$dates)],
                                         age = host_data_list[[i]]$demog_init$age)
  )
}
df_demog_kids <- df_demog %>% filter(age <= age_threshold)
df_demog_ad <- df_demog %>% filter(age > age_threshold)
n_kids = sum(age_vec <= age_threshold)
n_ad = sum(age_vec > age_threshold)
n_sim = max(titer_df_all$sim)

## Make plots -------------------------------------------------------------------------------------------------------------------------

## Distribution of n-fold titer rises:
fold_rises_plot_filename <- paste0("./Output_plots/fold_rises_",this_subtype,".pdf")
source("./Plotting_scripts/fold_rises_8_2018.R")
if(save_plots == T){
  save_plot(fold_rises_plot_filename, p_fold_rises , base_width = 6, base_height = 3)
}

## Coefficient of titer variation
CV_children_plot_filename <- paste0("./Output_plots/CV_children_",this_subtype,".pdf")
CV_adults_plot_filename <- paste0("./Output_plots/CV_adults_",this_subtype,".pdf")
source("./Plotting_scripts/CV_plot.R")
if(save_plots == T){
  save_plot(CV_children_plot_filename, p_fold_rises_children , base_width = 6, base_height = 3)
  save_plot(CV_adults_plot_filename, p_fold_rises_adults , base_width = 6, base_height = 3)
}

