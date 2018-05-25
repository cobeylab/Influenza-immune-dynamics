########################################################################
require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(lubridate)
require(dplyr)
source("../Utility_scripts/model_simulation_functions_multi_strain.R")
textSize = 12
source("../Utility_scripts/plot_themes.R")
require(pomp)
require(panelPomp)
require(cowplot)
select <- dplyr::select
summarize <- dplyr::summarize
n_strains = 1

## Input specs ## ----------------------------------------------------------------------------------------------------------------------------
model_versions = c("H1_kids","H1_adults","H3_kids","H3_adults")
sim_filenames = c("test_1000_sim_H1_kids.rda", "test_1000_sim_H1_adults.rda","test_1000_sim_H3_kids.rda","test_1000_sim_H3_adults.rda")
save_plots = F
## ----------------------------------------------------------------------------------------------------------------------------

make_boosting_plot <- function(titer_df_all){
      df_sim <- titer_df_all %>% 
        filter( date > date_inf) %>% 
        mutate(boost_obs_sim = floor(boost_obs_sim)) %>% 
        filter(boost_obs_sim < 10 & boost_obs_sim >=0) %>% 
        group_by(sim,boost_obs_sim) %>% 
        summarize(count = n()) %>% 
        group_by(boost_obs_sim) %>% 
        summarize(mean = mean(count),
                  LCI = quantile(count, c(.05,.95))[1],
                  UCI = quantile(count, c(.05,.95))[2]) %>% 
        mutate(data = "Simulated") %>% 
        rename(boost = boost_obs_sim)
        
      df_obs <- titer_df_all %>%  
        filter(date > date_inf & boost_obs >= 0) %>% 
        group_by(sim,boost_obs) %>% 
        summarize(count = n()) %>% 
        group_by(boost_obs) %>% 
        summarize(mean = median(count),
                  LCI = quantile(count, c(.05,.95))[1],
                  UCI = quantile(count, c(.05,.95))[2]) %>% 
        mutate(data = "Observed") %>% 
        rename(boost = boost_obs)
      
      
      
      df <- rbind(df_sim,df_obs)
      df_NA <- data.frame(boost = rep(c(0:10),2),
                          mean = NA, 
                          LCI = NA,
                          UCI = NA,
                          data = c(rep("Simulated", 11),rep("Observed", 11))
      ) 
      
      df <- rbind(df, df_NA)
      
      p <- ggplot(df, aes(x= boost, y = mean, fill = data)) + 
        geom_bar(stat = "identity",position = position_dodge(), width = 1) +
        geom_errorbar(aes(ymax= UCI, ymin = LCI), position = position_dodge(.9)) + 
        xlab("") +
        ylab("") + 
        plot_themes + 
        theme(legend.position = "none")
  return(p)
}

p_list <- list()
for(i in c(1:length(sim_filenames))){
  load(sim_filenames[i])
  p <- make_boosting_plot(titer_df_all)
  p_list[[i]] <- p
  rm(titer_df_all)
}

p_combined <- plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], ncol = 2, align = "hv")
if(save_plots){
  save_plot("./Output_plots/boosting_fig.pdf", p_combined, base_width = 183, base_height = 120, units = "mm")
}
