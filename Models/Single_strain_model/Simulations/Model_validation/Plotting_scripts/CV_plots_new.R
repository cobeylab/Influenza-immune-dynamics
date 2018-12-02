## Scripts to analyze the raw data 
library(RSQLite)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(corrplot)
library(cowplot)
library(lubridate)
textSize <- 12
source(".././Utility_scripts/plot_themes.R")
source(".././Utility_scripts/data_formatting_functions.R")
save_plots = 1
summarize <- dplyr::summarize
rename <- dplyr::rename
include_pilot = F
dbFilename <- ".././Data/hongkongflu.sqlite"
contains = dplyr::contains
select = dplyr::select

## Simulated data
titer_df_all_full <- titer_df_all
titer_df_all_full$h_obs = titer_df_all_full$h_obs + 1 
titer_df_all_full$h_obs_sim = titer_df_all_full$h_obs_sim + 1

make_plot <- function(titer_df, filename){
  test <- titer_df %>% 
    group_by(sim, ind) %>% 
    summarize(h_baseline_obs = as.numeric(h_obs[1] == 1),
              h_baseline_sim = as.numeric(h_obs_sim[1] == 1),
              #n_undetectable_obs = sum(as.numeric(h_obs[1] == 1)),
              #n_undetectable_sim = sum(as.numeric(h_obs[1] == 1)),
              #n_detectable_obs = sum(as.numeric(h_obs[1] > 1)),
              cv_obs = sd(h_obs)/mean(h_obs),
              cv_sim = sd(h_obs_sim)/mean(h_obs_sim)) #%>% 
    #group_by(ind) %>% 
    #summarize(mean_cv_sim = mean(cv_sim),
           #   sd_cv_sim = sd(cv_sim),
           #   mean_cv_obs = mean(cv_obs),
           #   sd_cv_obs = 0)
  
  test_obs <- test %>% 
    #filter(sim ==1 ) %>% 
    select(c(ind, contains("obs"))) %>%
    mutate(data = "Observed") %>% 
    rename(mean = cv_obs, 
           baseline = h_baseline_obs)
                   
  test_sim <- test %>% 
    select(c(ind, contains("sim"))) %>%
    mutate(data = "Simulated") %>% 
    rename(mean = cv_sim, 
           baseline = h_baseline_sim)
  
  df = rbind(test_obs, test_sim)
  
  df2 <- df %>% 
    group_by(data,sim,baseline) %>% 
    mutate(mean = cut(mean,30)) %>% 
    ungroup %>% 
    group_by(data,sim,baseline,mean) %>% 
    summarize(n = length(ind)) %>% 
    ungroup %>% 
    as.data.frame %>% 
    group_by(data,sim,baseline) %>% 
    mutate(n_ind = sum(n)) %>% 
    ungroup %>% 
    as.data.frame %>% 
    mutate(frac = n/n_ind) %>% 
    group_by(data,sim,baseline,mean) %>% 
    summarize(mean_frac = mean(frac))
  
  
  df$baseline_f = factor(df$baseline, labels = c("Detectable initial titer", "Undetectable intiial titer"))
  
  # df1 <- df %>% filter(data == "Observed" & baseline == 1)
  # df2 <- df %>% filter(data == "Observed" & baseline == 0)
  # df3 <- df %>% filter(data == "Simulated" & baseline == 1)
  # df4 <- df %>% filter(data == "Simulated" & baseline == 0)
  # p <- ggplot(df, aes(x= mean)) + 
  #   #geom_bar(stat = "identity", position = "dodge", aes(color = data), alpha = .4) + 
  #   #stat_bin(aes(y=..density.., fill = data),alpha = .4, position='dodge')+
  #   geom_histogram(aes(y = ..count../sum(..count..), fill = data), alpha = .4, position = "dodge") + 
  #   geom_histogram(data = df1, aes(y = ..count../sum(..count..), fill = data), alpha = .4, position = "dodge") + 
  #   geom_histogram(data = df2, aes(y = ..count../sum(..count..), fill = data), alpha = .4, position = "dodge") + 
  #   geom_histogram(data = df3, aes(y = ..count../sum(..count..), fill = data), alpha = .4, position = "dodge") + 
  #   facet_wrap(~ baseline_f) + 
  #   xlab("") + 
  #   ylab("") + 
  #   plot_themes 
  
  
  p <- ggplot(df, aes(x= mean)) + 
    #geom_bar(stat = "identity", position = "dodge", aes(color = data), alpha = .4) + 
    stat_bin(aes(y=..density.., fill = data),alpha = .4, position='dodge') + 
    xlab("") + 
    ylab("") + 
    facet_wrap(~baseline_f) + 
    plot_themes 


  save_plot(filename, p , base_width = 6, base_height = 3)
  
}

## Plot for adults
make_plot(titer_df = titer_df_all_full %>% filter(age > 20), filename = adults_CV_plot_filename)

## Plot for kids
make_plot(titer_df = titer_df_all_full %>% filter(age < 20), filename = kids_CV_plot_filename)

