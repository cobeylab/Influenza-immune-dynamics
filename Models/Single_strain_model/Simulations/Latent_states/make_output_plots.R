########################################################################
require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(lubridate)
require(dplyr)
source("../../Utility_scripts/model_functions.R")
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

# Load in the simulated host list and simulated output for adults and children 
pomp_object_filename_adults <- "" # Name of file (".rda") containing pomp object for adults with extra simulated observations
pomp_object_filename_kids <- "" # Name of file (".rda") containing pomp object for children with extra simulated observations
model_predictions_filename_adults <- "" # Name of file (".rda") containing model simulations for adults
model_predictions_filename_kids <- "" # Name of file (".rda") containing model simulations for adults

load(pomp_object_filename_adults)
rm(panelObject)
load(model_predictions_filename_adults)
titers_adults <- titer_df_all
titers_adults$child = 0
infections_adults <- inf_times_df_all
infections_adults$child = 0
hosts_adults <- host_data_list

load(pomp_object_filename_kids)
rm(panelObject)
load(model_predictions_filename_kids)
titers_children <- titer_df_all
titers_children$child = 1
infections_children <- inf_times_df_all
infections_children$child = 1
hosts_children <- host_data_list

inf_times_df_all <- rbind(infections_children, infections_adults)
titer_df_all <- rbind(titers_children, titers_adults)
host_data_list <- c(hosts_children, hosts_adults)


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
for( i in c(1:length(host_data_list))){
  age_vec[i] <- as.numeric( host_data_list[[i]]$demog_init["age"])
  df_demog <- rbind(df_demog, data.frame(ind = i,
                                         first_date = host_data_list[[i]]$dates[1],
                                         last_date = host_data_list[[i]]$dates[length(host_data_list[[i]]$dates)],
                                         age = host_data_list[[i]]$demog_init["age"])
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

if(save_plots){
  save_plot(latent_susceptibility_plot_filename, p_latent_susceptibility, base_width = 6, base_height = 4 )
}
