## Plot distribution of fold rises ## ---------------------------------------------------------------------

# N-fold rises 
#rm(df_sim,df_obs,df2,df4,df8,df)

# 2 - fold 
    df_sim <- titer_df_all %>% 
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
                min = min(n),
                max = max(n),
                data = "sim") %>% 
      rename(n = n_2)
    
    df_obs <- titer_df_all %>% 
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
                min = min(n),
                max = max(n),
                data = "obs") %>% 
      rename(n=n_2)
      
    df2 <- rbind(df_sim,df_obs) %>% 
      mutate(n_fold = "2 - fold")
# 4-fold
    df_sim <- titer_df_all %>% 
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
                min = min(n),
                max = max(n),
                data = "sim") %>% 
      rename(n = n_4)
    
    df_obs <- titer_df_all %>% 
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
                min = min(n),
                max = max(n),
                data = "obs") %>% 
      rename(n = n_4)
    
    df4 <- rbind(df_sim,df_obs) %>% 
      mutate(n_fold = "4 - fold")
## 8-fold 
    df_sim <- titer_df_all %>% 
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
    
    df_obs <- titer_df_all %>% 
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

p = ggplot(df,aes(x=n, y = mean)) + 
  geom_point(aes(color =data)) + 
  geom_line(aes(x=n,y=min,color = data)) + 
  geom_line(aes(x=n,y=max,color=data)) + 
  facet_wrap(~n_fold)

n_ind = length(unique(titer_df_all$ind))
p_fold_rises <-  ggplot(df %>% filter(data == "sim"), aes( x = n , y = mean/n_ind)) + 
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
  facet_wrap(~n_fold, scales = "free_x") + 
  plot_themes 

