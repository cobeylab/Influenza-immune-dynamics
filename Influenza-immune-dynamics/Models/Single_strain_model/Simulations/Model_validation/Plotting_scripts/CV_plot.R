##-----------------------------------------------------------------------------------------------------
## Sylvia Ranjeva, May 2018 
## Plot the coefficient of titer variation among individuals 
## -------------------------------------------------------------------------------------------------------------------
make_cv_plot <- function(titer_df_sub){
  cv_data <- titer_df_sub %>% 
    group_by(sim, ind) %>% 
    summarize(h_baseline_obs = as.numeric(h_obs[1] == 1),
              h_baseline_sim = as.numeric(h_obs_sim[1] == 1),
              cv_obs = sd(h_obs)/mean(h_obs),
              cv_sim = sd(h_obs_sim)/mean(h_obs_sim)) 
  
  cv_data_obs <- cv_data %>% 
    select(c(ind, contains("obs"))) %>%
    mutate(data = "Observed") %>% 
    rename(mean = cv_obs, 
           baseline = h_baseline_obs)
                   
  cv_data_sim <- cv_data %>% 
    select(c(ind, contains("sim"))) %>%
    mutate(data = "Simulated") %>% 
    rename(mean = cv_sim, 
           baseline = h_baseline_sim)
  
  df = rbind(cv_data_obs, cv_data_sim)
  df$baseline_f = factor(df$baseline, labels = c("Detectable initial titer", "Undetectable intiial titer"))
  
  p_cv <- ggplot(df, aes(x= mean)) + stat_bin(aes(y=..density.., fill = data),alpha = .4, position='dodge')+
    facet_wrap(~ baseline_f) + 
    xlab("CV") + 
    ylab("Density") + 
    plot_themes 
  
  return(p_cv)
}

p_cv_children <- make_cv_plot(titer_df_sub = titer_df_all %>% filter(child == 1))
p_cv_adults <- make_cv_plot(titer_df_sub = titer_df_all %>% filter(child == 0))