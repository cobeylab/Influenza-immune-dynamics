contact_matrix <- c(3.11, 0.51, 0.14, 0.07, 0.67, 0.86, 1.58, 0.85, 0.71, 0.21, 0.34, 0.50, 0.46,
              0.13, 0.05, 0.43, 3.50, 0.77, 0.27, 0.22, 0.48, 0.86, 0.90, 0.98, 0.48, 0.15,
              0.06, 0.22, 0.17, 0.27, 0.05, 0.42, 6.56, 0.37, 0.47, 0.28, 0.61, 0.72, 1.18,
              0.73, 0.61, 0.06, 0.11, 0.07, 0.37, 0.02, 0.05, 0.49, 2.35, 1.47, 0.23, 0.69,
              0.23, 0.63, 0.71, 0.79, 0.16, 0.05, 0.00, 0.04, 0.08, 0.07, 0.16, 0.53, 2.16,
              0.95, 1.54, 0.68, 1.37, 0.45, 1.04, 0.89, 0.49, 0.02, 0.00, 0.16, 0.02, 0.01,
              0.59, 0.45, 1.43, 0.76, 0.38, 0.65, 0.32, 0.17, 0.21, 0.40, 0.23, 0.16, 0.38,
              0.28, 0.84, 0.18, 0.80, 1.03, 2.32, 0.58, 0.92, 0.49, 0.51, 0.34, 0.79, 0.12,
              0.23, 0.45, 0.31, 0.23, 0.14, 0.23, 0.40, 0.78, 0.86, 1.08, 0.33, 0.33, 0.17,
              0.28, 0.13, 0.17, 0.12, 0.23, 0.42, 0.13, 0.32, 0.64, 0.76, 0.65, 1.15, 0.86, 
              0.70, 0.31, 0.22, 0.05, 0.08, 0.06, 0.19, 0.30, 0.38, 0.41, 0.62, 0.99, 0.84,
              1.14, 0.98, 0.82, 0.50, 0.53, 0.01, 0.14, 0.07, 0.14, 0.26, 0.34, 0.49, 0.55,
              0.84, 0.47, 0.77, 1.25, 0.97, 0.51, 0.33, 0.07, 0.20, 0.47, 0.06, 0.08, 0.16,
              0.41, 0.51, 0.66, 0.63, 0.82, 0.70, 1.18, 0.73, 0.65, 0.41, 0.31, 0.10, 0.03,
              0.04, 0.04, 0.23, 0.36, 0.40, 0.38, 0.30, 0.42, 0.51, 0.80, 0.61, 0.35, 0.18,
              0.21, 0.00, 0.09, 0.05, 0.23, 0.34, 0.41, 0.51, 0.57, 0.19, 0.32, 0.39, 0.67,
              0.26, 0.18, 0.03, 0.02, 0.00, 0.22, 0.16, 0.31, 0.38, 0.54, 0.54, 0.26, 0.19,
              0.33, 0.55, 0.44, 0.66)

age_vec <- seq(0,70, by = 5)
contact_df <- cbind(expand.grid(age_index = age_vec, age_contact = age_vec), n_contacts = contact_matrix)
mapping_vec <- c(0,10,20,40,65)
#contact_df$group_index <- NA
contact_df$group_contact <- NA
contacts_index_total <- data.frame(age_participant = c(0,10,20,40,65),
                      total = c(1,.953,.746,.751,.562)*7.65,
                      percent_ILI = c(.19, .33, .30, .14, .04))

ILI_dist_vec <- data.frame(ages = c(0,10,20,40,65), percent = c(.19, .33, .30, .14, .04))

for(i in c(1:nrow(contact_df))){
  #mapping_ind <- min(which(contact_df[i,]$age_index <= mapping_vec))
  mapping_contact <- max(which(contact_df[i,]$age_contact >= mapping_vec))
  #contact_df[i,]$group_index <- mapping_vec[mapping_ind]
  contact_df[i,]$group_contact <- mapping_vec[mapping_contact]
}
contact_df_reduced <- contact_df %>% 
  group_by(age_index, group_contact) %>% 
  summarize(mean_contacts = mean(n_contacts)) %>% 
  as.data.frame() %>% 
  group_by(age_index) %>% 
  mutate(frac_contacts = mean_contacts/sum(mean_contacts))

contact_df_reduced$daily_contacts_index = NA
contact_df_reduced$ILI_weight = NA
for(i in c(1:nrow(contact_df_reduced))){
  ind = max(which(contacts_index_total$age_participant <= contact_df_reduced$age_index[i]))
  ind_group = which(contact_df_reduced$group_contact[i] == ILI_dist_vec$ages)
  contact_df_reduced$daily_contacts_index[i] = contacts_index_total[ind,]$total
  contact_df_reduced$ILI_weight[i] = ILI_dist_vec[ind_group,]$percent
}

contact_df_red2 <- contact_df_reduced %>% 
  mutate(weighted_daily_contacts = daily_contacts_index*frac_contacts) %>% 
  #select(age_index,group_contact, daily_contacts_index, weighted_daily_contacts,ILI_weight) %>% 
  mutate(frac_weighted_daily_contacts = weighted_daily_contacts/sum(weighted_daily_contacts)) %>% 
  mutate(ILI_weighted_contact_frac = frac_weighted_daily_contacts*ILI_weight)
  
contacts <- contact_df_red2 %>% 
  group_by(age_index) %>% 
  summarize(tot_weighted_contacts = sum(ILI_weighted_contact_frac)*length(ILI_weighted_contact_frac)*unique(daily_contacts_index)) %>% 
  as.data.frame() %>% 
  rename(age_participant = age_index,
         total = tot_weighted_contacts)
