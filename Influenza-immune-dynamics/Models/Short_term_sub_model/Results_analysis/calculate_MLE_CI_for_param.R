
### ---------------------------------------------------------------------------------------
## Sylvia Ranjeva, May 2018
## From the profile for one parameter, calculate the MLE and the 95% CI
## ---------------------------------------------------------------------------------------
require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../../Utility_scripts/model_functions.R")
textSize = 12
source("../../Utility_scripts/plot_themes.R")
require(pomp)
require(panelPomp)
select <- dplyr::select
summarize <- dplyr::summarise 

## Input specs ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
profile_filename <- "" # Name of the file (".csv") containing the output for the profile likelihood search
param_name = "" # Name of the focal parameter
n_bin_points = 200 # Number of points over which to profile
alpha_min =  # Minimum focal parameter value
alpha_min =  # Maximum focal parameter value
bin_points <- seq(alpha_min, alpha_max, length.out = n_bin_points)

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("param_names.rda")
results <- as.data.frame(read.csv( profile_filenames[i] ))%>% select(-se) 
names(results) <- c(names, "loglik", "loglik_se", "n_mif", "n_part", "chain")

results_profiles <- data.frame()
this_param <- param_name
results$focal_param <- results[,names(results) == this_param]
results_sub <- results %>% group_by(focal_param) %>%  summarize(LL = max(loglik)) %>% mutate(dLL = LL - max(LL))
results_sub$param = this_param
results_profiles <- rbind(results_profiles,results_sub)
results_sub$dLL <- results_sub$LL - max(results_sub$LL)

results_profiles = data.frame()
res = results
df_binned <- data.frame()
for( j in c(1:(length(bin_points)-1))){
  df_sub <- res %>% filter(focal_param >= bin_points[j] & focal_param < bin_points[j+1]) %>% 
    arrange(-loglik)
  if(nrow(df_sub) > 0){
    df_binned <- rbind(df_binned,df_sub[1,])
  }
}
df_binned$param = param_name
results_profiles = rbind(results_profiles,df_binned)
results_profiles$dLL = results_profiles$loglik - max(results_profiles$loglik)


df_MLE <- data.frame()

sub <- results_profiles 
mcap_sub <- mcap(lp = sub$dLL, param = sub$focal_param, lambda = .75)
fit = mcap_sub$fit
df <- data.frame(MLE = mcap_sub$mle,
                   LCI = mcap_sub$ci[1],
                   UCI = mcap_sub$ci[2],
                   param = param_name)
