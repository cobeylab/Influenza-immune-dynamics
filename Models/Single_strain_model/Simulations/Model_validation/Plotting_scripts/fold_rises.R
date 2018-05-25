
##-----------------------------------------------------------------------------------------------------
## Sylvia Ranjeva, May 2018 
## Plot the distribtion of n-fold titer rises among individuals
## ----------------------------------------------------------------------------------------------------------------

## Function to apply measurement error to latent titer ## ---------------------------------------------------------------
err_func = function(x,thres,sig,sig2){
  if(x > thres ){
    y = floor(rnorm(1,mean = x,sd = sig))
    if(y > 10){
      y = 10
    }
    if(y < 2){
      y = 1
    }
  }
  if(x <= thres){
      sigma = sig2
    y = floor(rnorm(1,mean = x,sd = sigma))
    if(y > 10){
      y = 10
    }
    if(y < 2){
      y = 1
    }
  }
  return(y)
}
titer_df_all$h_obs_sim = sapply(titer_df_all$h_latent, err_func2, thres = 1, sig = 1.29, sig2 = 0.74) # Apply the measurement error

## Calculate the distribution of 2-fold, 4-fold, and 8-fold titer rises ## ---------------------------------------------------------------
make_fold_rises_plot <- function(titer_df_sub){
  #rm(df,df2,df4,df8,df_sim,df_obs) # Make sure these data frames do not already exist (I've made this mistake more than once...)

  # 2 - fold 
    df_sim <- titer_df_sub %>% 
      select(ind,sim,time,date,age,h_obs_sim) %>% 
      group_by(ind,sim,age) %>% 
      mutate(diff = c(NA,diff(h_obs_sim)),
             max = max(h_obs_sim)) %>% 
      filter(!is.na(diff)) %>% 
      mutate(twofold = as.numeric(diff > 0),
             fourfold = as.numeric(diff > 1),
             eightfold = as.numeric(diff > 2)) %>% 
      group_by(ind,sim,age,max) %>% 
      summarize(n_2 = sum(twofold),
                n_4 = sum(fourfold),
                n_8 = sum(eightfold)) %>% 
      group_by(sim,n_2) %>% 
      summarize(n =length(ind)) %>% 
      group_by(n_2) %>% 
      summarize(mean = mean(n),
                min = min(n),#quantile(n,c(.025,.975))[1],#min(n),
                max = max(n),
                data = "sim") %>% 
      rename(n = n_2)
    
    df_obs <- titer_df_sub %>% 
      filter(sim ==1) %>% 
      select(ind,sim,time,date,age,h_obs) %>% 
      group_by(ind,sim,age) %>% 
      mutate(diff = c(NA,diff(h_obs))) %>% 
      filter(!is.na(diff)) %>% 
      mutate(twofold = as.numeric(diff > 0),
             fourfold = as.numeric(diff > 1),
             eightfold = as.numeric(diff > 2)) %>% 
      group_by(ind,sim,age) %>% 
      summarize(n_2 = sum(twofold),
                n_4 = sum(fourfold),
                n_8 = sum(eightfold)) %>% 
      group_by(sim,n_2) %>% 
      summarize(n =length(ind)) %>% 
      group_by(n_2) %>% 
      summarize(mean = mean(n),
                min = min(n),#quantile(n,c(.1,.9))[1],#min(n),
                max = max(n), #quantile(n,c(.1,.9))[2],
                data = "obs") %>% 
      rename(n=n_2)
      
    df2 <- rbind(df_sim,df_obs) %>% 
      mutate(n_fold = "2 - fold")
  # 4-fold
    df_sim <- titer_df_sub %>% 
      select(ind,sim,time,date,age,h_obs_sim) %>% 
      group_by(ind,sim,age) %>% 
      mutate(diff = c(NA,diff(h_obs_sim))) %>% 
      filter(!is.na(diff)) %>% 
      mutate(twofold = as.numeric(diff > 0),
             fourfold = as.numeric(diff > 1),
             eightfold = as.numeric(diff > 2)) %>% 
      group_by(ind,sim,age) %>% 
      summarize(n_2 = sum(twofold),
                n_4 = sum(fourfold),
                n_8 = sum(eightfold)) %>% 
      group_by(sim,n_4) %>% 
      summarize(n =length(ind)) %>% 
      group_by(n_4) %>% 
      summarize(mean = mean(n),
                min = min(n),#quantile(n,c(.025,.975))[1],#min(n),
                max = max(n),#quantile(n,c(.025,.975))[2],
                data = "sim") %>% 
      rename(n = n_4)
    
    df_obs <- titer_df_sub %>% 
      filter(sim == 1) %>% 
      select(ind,sim,time,date,age,h_obs) %>% 
      group_by(ind,sim,age) %>% 
      mutate(diff = c(NA,diff(h_obs))) %>% 
      filter(!is.na(diff)) %>% 
      mutate(twofold = as.numeric(diff > 0),
             fourfold = as.numeric(diff > 1),
             eightfold = as.numeric(diff > 2)) %>% 
      group_by(ind,sim,age) %>% 
      summarize(n_2 = sum(twofold),
                n_4 = sum(fourfold),
                n_8 = sum(eightfold)) %>% 
      group_by(sim,n_4) %>% 
      summarize(n =length(ind)) %>% 
      group_by(n_4) %>% 
      summarize(mean = mean(n),
                min = min(n),#quantile(n,c(.025,.975))[1],#min(n),
                max = max(n),#quantile(n,c(.025,.975))[2],
                data = "obs") %>% 
      rename(n = n_4)
    
    df4 <- rbind(df_sim,df_obs) %>% 
      mutate(n_fold = "4 - fold")
  # 8-fold 
    df_sim <- titer_df_sub %>% 
      select(ind,sim,time,date,age,h_obs_sim) %>% 
      group_by(ind,sim,age) %>% 
      mutate(diff = c(NA,diff(h_obs_sim))) %>% 
      filter(!is.na(diff)) %>% 
      mutate(twofold = as.numeric(diff > 0),
             fourfold = as.numeric(diff > 1),
             eightfold = as.numeric(diff > 2)) %>% 
      group_by(ind,sim,age) %>% 
      summarize(n_2 = sum(twofold),
                n_4 = sum(fourfold),
                n_8 = sum(eightfold)) %>% 
      group_by(sim,n_8) %>% 
      summarize(n =length(ind)) %>% 
      group_by(n_8) %>% 
      summarize(mean = mean(n),
                min = min(n),
                max = max(n),
                data = "sim") %>% 
      rename(n = n_8)
    
    df_obs <- titer_df_sub %>% 
      filter(sim == 1) %>% 
      select(ind,sim,time,date,age,h_obs) %>% 
      group_by(ind,sim,age) %>% 
      mutate(diff = c(NA,diff(h_obs))) %>% 
      filter(!is.na(diff)) %>% 
      mutate(twofold = as.numeric(diff > 0),
             fourfold = as.numeric(diff > 1),
             eightfold = as.numeric(diff > 2)) %>% 
      group_by(ind,sim,age) %>% 
      summarize(n_2 = sum(twofold),
                n_4 = sum(fourfold),
                n_8 = sum(eightfold)) %>% 
      group_by(sim,n_8) %>% 
      summarize(n =length(ind)) %>% 
      group_by(n_8) %>% 
      summarize(mean = mean(n),
                min = min(n),
                max = max(n),
                data = "obs") %>% 
      rename(n = n_8)
    
    df8 <- rbind(df_sim,df_obs) %>% 
      mutate(n_fold = "8 - fold")
    
  df = rbind(df2,df4,df8)    

  n_ind = length(unique(titer_df_sub$ind))
  p_frac <-  ggplot(df %>% filter(data == "sim"), aes( x = n , y = mean/n_ind)) + 
    geom_point(aes(color = data)) + 
    geom_line(aes(color = data), linetype =2 ) + 
    geom_ribbon(aes(ymin = min/n_ind, ymax = max/n_ind), 
                fill = "blue",
                linetype = 2,
                alpha=.15) + 
    geom_point(data = df %>% filter(data == "obs"), aes( x = n, y = mean/n_ind, color = data)) +
    geom_line(data = df %>% filter(data == "obs"), aes( x = n, y = mean/n_ind, color = data), linetype =2) +
    xlab("Number of n-fold rises over followup") + 
    labs(color = "Data type") + 
    ylab("Fraction of population") +
    facet_wrap(~n_fold) + 
    xlim(0,4) + # No individuals with > 4 observed n-fold rises 
    plot_themes 

  return(p_frac)
}

p_fold_rises_adults <- make_fold_rises_plot(titer_df_sub = titer_df_all %>% filter(child == 0))
p_fold_rises_children <- make_fold_rises_plot(titer_df_sub = titer_df_all %>% filter(child == 1))
