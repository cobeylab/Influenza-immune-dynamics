simulation_timestep = 2.5 #Days
dat_all_sim = data.frame()

for( i in c(1:n_sim)){
    inf <- inf_times_df_all %>% 
      filter(sim==i) %>% 
      group_by(ind, child) %>% 
      mutate(inf_num = c(1:length(inf_time)),
      max_inf = length(inf_time)) %>% 
      as.data.frame %>% 
      select(-c(age,date,sim))
  
    dat <- titer_df_all %>% 
      filter(sim ==i ) %>%
      left_join(.,inf, by = c("ind","child")) %>% 
      mutate(dt = time - inf_time) %>% 
      filter(max_inf ==1 & dt > 0)

  dat_all_sim = rbind(dat_all_sim,dat)
}

dat_all_sim$child_f <- factor(dat_all_sim$child, labels = c("Adults","Children"))
p_latent_susceptibility <- ggplot() + 
  geom_line(data = dat_all_sim, aes(x= dt/365, y = q_latent, group = interaction(ind,sim)),linetype = 1, alpha = .1) + 
  plot_themes + 
  xlab("Time since infection (y)") + ylab("Susceptibility") + 
  xlim(2*simulation_timestep/365,4) + 
  ylim(0,1) + 
  facet_wrap(~child_f)



