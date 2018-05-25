
df_all <- data.frame()
 
for(j in c(1:n_sim)){
  cat("j is ", j, "\n")
  test_sim <- inf_times_df_all %>% filter(sim == j)
  test_titer <- titer_df_all %>% filter(sim == j)
  test_sim$child <- as.factor(test_sim$age < 15)
  test_titer$child <- as.factor(test_titer$age < 15)
  test_sim_kids = test_sim %>% filter(child ==T)
  test_sim_ad = test_sim %>% filter(child ==F)
  test_titer_kids = test_titer %>% filter(child ==T)
  test_titer_ad = test_titer %>% filter(child ==F)
  AR_kids <- array(NA, length(test_seq))
  AR_ad <- array(NA,length(test_seq))
  for(i in c(1:(length(test_seq)-1))){
     n_kids <- sum(as.numeric(as.Date(df_demog_kids$first_date,origin = '1970-1-1')) < as.numeric(as.Date(test_seq[i],origin = '1970-1-1')) & as.numeric(as.Date(df_demog_kids$last_date,origin = '1970-1-1')) > as.numeric(as.Date(test_seq[i+1],origin = '1970-1-1')) ,na.rm=T)
    n_ad <- sum(as.numeric(as.Date(df_demog_ad$first_date,origin = '1970-1-1')) < as.numeric(as.Date(test_seq[i],origin = '1970-1-1')) & as.numeric(as.Date(df_demog_ad$last_date,origin = '1970-1-1')) > as.numeric(as.Date(test_seq[i+1],origin = '1970-1-1')) ,na.rm=T)
    AR_kids[i] <- sum(test_sim_kids$date > test_seq[i] & test_sim_kids$date < test_seq[i+1])/n_kids
    AR_ad[i] <- sum(test_sim_ad$date > test_seq[i] & test_sim_ad$date < test_seq[i+1])/n_ad
  }
  df <- rbind(data.frame(date = test_seq, AR = AR_kids,  age = "Child"), 
              data.frame(date = test_seq, AR = AR_ad,  age = "Adult"))

  df$sim = j
  df_all <- rbind(df_all,df)
}


df_AR <- df_all %>% filter(!is.na(AR)) %>% group_by(date, age) %>% summarize(mean = mean(AR), 
                                                                        LCI = quantile(AR,c(.05,.95))[1],
                                                                        UCI = quantile(AR,c(.05,.95))[2]
                                                                        )

dat_sub <- intensity_data %>% filter(full_date >=min(df$date) & full_date <= max(df$date))
dat_sub$month = month(dat_sub$full_date)
dat_sub$year = year(dat_sub$full_date)
dat_sub_monthly <- dat_sub %>% group_by(month,year) %>% summarize(L_monthly = max(sum(L),1e-5))
dat_sub_monthly$time = as.Date(paste0(dat_sub_monthly$year,"-",dat_sub_monthly$month, "-", 1), origin = "1970-1-1")


p_community_intensity <- ggplot() + 
  geom_line(data = dat_sub_monthly, aes(x = time, y = L_monthly), linetype = 2) + 
  xlab("") + 
  plot_themes + 
  ylab("Monthly intensity") 

p_AR <- ggplot(df_AR, aes(x = date, y = mean)) + geom_line(aes(color=age)) + 
  geom_ribbon(aes(ymin=LCI,ymax=UCI, fill = age), alpha="0.25") +
  ylab("Monthly incidence") + 
  xlab("") + 
  plot_themes + 
  theme(legend.position = c(0.9, 0.75))

p <- plot_grid(p_community_intensity, p_AR, ncol = 1, rel_heights = c(1,2), align = "hv")


