########################################################################
## From the results of the model, analyze attack rates etc
## Attack rate, population dist. of ab by age at specified time points 
########################################################################
require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
require(pomp)
require(panelPomp)
require(cowplot)
source("../../Utility_scripts/model_simulation_functions_multi_strain.R")
textSize = 12
source("../../Utility_scripts/plot_themes.R")
select <- dplyr::select
n_strains = 1

if(this_subtype == "H3N2"){
  #H3N2 epidemics 
  ep1start <- as.Date("2010-8-1",origin = "1970-1-1")
  ep1end <- as.Date("2010-10-16",origin = "1970-1-1")
  ep2start <-as.Date("2012-3-11",origin = "1970-1-1")
  ep2end <-  as.Date("2012-6-23",origin = "1970-1-1")
  ep3start <- as.Date("2013-7-14",origin = "1970-1-1")
  ep3end <- as.Date("2013-10-12",origin = "1970-1-1")
}

if(this_subtype == "pH1N1"){
  #PH1N1 epidemics
  ep1start <- as.Date("2009-11-1",origin = "1970-1-1")
  ep1end <- as.Date("2010-4-16",origin = "1970-1-1")
  ep2start <- as.Date("2010-12-9",origin = "1970-1-1")
  ep2end <- as.Date("2011-2-26",origin = "1970-1-1")
  ep3start <- as.Date("2013-1-1",origin = "1970-1-1")
  ep3end <- as.Date("2013-4-20",origin = "1970-1-1")
}

sampled_ep1 <- df_demog %>% filter(first_date<ep1start & last_date > ep1end) %>% select(ind) %>% unlist %>% as.numeric
sampled_ep2 <- df_demog %>% filter(first_date<ep2start & last_date > ep2end) %>% select(ind) %>% unlist %>% as.numeric
sampled_ep3 <- df_demog %>% filter(first_date<ep3start & last_date > ep3end) %>% select(ind) %>% unlist %>% as.numeric

n_sim <- length(unique(inf_times_df_all$sim))
df_all <- data.frame()
df_lm = data.frame()
inf_times_df_all$child <- as.factor(inf_times_df_all$age <=15)
for( j in c(1:n_sim)){
  dat <- inf_times_df_all %>% filter(sim==j)
  ep1_dat <- dat %>% filter(date >= ep1start & date <= ep1end) %>% 
    select(ind) %>% 
    unlist
  ep2_dat <- dat %>% filter(date >= ep2start & date <= ep2end) %>% 
    select(ind) %>% 
    unlist
  ep3_dat <- dat %>% filter(date >= ep3start & date <= ep3end) %>% 
    select(ind) %>% 
    unlist
  df <- data.frame(ind = c(1:nrow(df_demog)), inf1 = 0, sampled_ep1 = 0, inf2 = 0 , sampled_ep2 = 0, inf3 = 0, sampled_ep3 = 0)
  for( i in c(1:nrow(df))){
    df[i,]$inf1 = as.numeric(i %in% ep1_dat)
    df[i,]$sampled_ep1 = as.numeric( i %in% sampled_ep1)
    df[i,]$inf2 = as.numeric(i %in% ep2_dat)
    df[i,]$sampled_ep2 = as.numeric( i %in% sampled_ep2)
    df[i,]$inf3 = as.numeric(i %in% ep3_dat)
    df[i,]$sampled_ep3 = as.numeric( i %in% sampled_ep3)
  }
  
  df_lm = rbind(df_lm, df)
  df_sim <- data.frame(inf1 = sum(df$inf1==1),
                       uninf1 = sum(df$inf1 == 0),
                       inf_1_uninf_2 = sum(df$inf1 == 1 & df$inf2 == 0 & df$sampled_ep2 ==1), 
                       inf_1_inf_2 = sum(df$inf1 == 1 & df$sampled_ep1 ==1 & df$inf2 == 1 & df$sampled_ep2 ==1),
                       uninf_1_uninf_2 = sum(df$inf1 == 0 & df$sampled_ep1 ==1 & df$inf2 == 0 & df$sampled_ep2 ==1),
                       uninf_1_inf_2 = sum(df$inf1 == 0 & df$sampled_ep1 ==1 & df$inf2 == 1),
                       inf_2_inf_3=  sum(df$inf1 == 1 & df$inf3 == 1),
                       inf_2_uninf_3 =  sum(df$inf2 == 1 & df$inf3 == 0 & df$sampled_ep3 ==1),
                       uninf_2_uninf_3 = sum(df$inf2 == 0 & df$sampled_ep2 ==1 & df$inf3 == 0 &  df$sampled_ep3 ==1),
                       uninf_2_inf_3 = sum(df$inf2 == 0 & df$sampled_ep2 ==1 & df$inf3 == 1),
                       uninf_1_uninf_3 = sum(df$inf1 == 0 & df$sampled_ep1 ==1 & df$inf3 == 0 & df$sampled_ep3 ==1),
                       inf_1_inf_3 =sum(df$inf1 == 1 & df$inf3 == 1),
                       inf_1_uninf_3 = sum(df$inf1 == 1 & df$inf3 == 0 & df$sampled_ep3 ==1),
                       uninf_1_inf_3 = sum(df$inf1 == 0 &  df$sampled_ep1 ==1 & df$inf3 == 1),
                       sim = j
  ) %>% 
    mutate( OR_1_2 = (inf_1_inf_2/uninf_1_inf_2)/(inf_1_uninf_2/uninf_1_uninf_2)) %>% 
    mutate(OR_2_3 = (inf_2_inf_3/uninf_2_inf_3)/(inf_2_uninf_3/uninf_2_uninf_3)) %>% 
    mutate(OR_1_3 = (inf_1_inf_3/uninf_1_inf_3)/(inf_1_uninf_3/uninf_1_uninf_3))
  df_all <- rbind(df_all,df_sim)
}


lm_1_2 = glm(inf2 ~ inf1 + (1|sim) + (1|ind),family= binomial(link='logit'),data=df_lm)
lm_2_3 = glm(inf3 ~ inf2 + (1|sim) + (1|ind),family= binomial(link='logit'),data=df_lm)
lm_1_3 = glm(inf3 ~ inf1 (1|sim) + (1|ind),family= binomial(link='logit'),data=df_lm)