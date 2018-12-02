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
pomp_object_filename <- "" # Name of file (".rda") containing pomp object with extra simulated observations
model_predictions_filename <- "" # Name of file (".rda") containing model simulations 

load(pomp_object_filename)
rm(panelObject)
load(model_predictions_filename)

titer_df_all$age = titer_df_all$age + titer_df_all$time/365
inf_times_df_all$child = as.numeric(inf_times_df_all$age < 15)
inf_times_df_all$adult = as.numeric(inf_times_df_all$age > 20)
titer_df_all$child = as.numeric(titer_df_all$age < 15)
titer_df_all$adult = as.numeric(titer_df_all$age > 20)

## Format data ## ----------------------------------------------------------------------------------------------------------------------------

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
test_seq <- seq.Date(as.Date("2009-12-1"), as.Date("2014-1-1"),by = "month")
n_sim = max(titer_df_all$sim)

## Make plots -------------------------------------------------------------------------------------------------------------------------

## Attack rate plo for this subtype:
attack_plot_filename <- paste0("./Output_plots/AR_plot_",this_subtype,".pdf")
source("./Plotting_scripts/attack_rate_plot.R")

if(save_plots){
  save_plot(attack_plot_filename, p, base_height = 6, base_width = 8)
}

## Plot of latent susceptibility after infection for adults and children for this subtype
latent_susceptibility_plot_filename <- paste0("./Output_plots/latent_susceptibility_",this_subtype,".pdf")
source("./Plotting_scripts/latent_susceptibility_after_infection_plot.R")

if(save_plots){
  save_plot(latent_susceptibility_plot_filename, p_latent_susceptibility, base_width = 6, base_height = 4 )
}
