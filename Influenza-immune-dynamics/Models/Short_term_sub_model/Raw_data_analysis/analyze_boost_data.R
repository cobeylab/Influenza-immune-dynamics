require(ggplot2)
require(cowplot)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../../../Utility_scripts/model_functions.R")
textSize = 12
source("../../../Utility_scripts/plot_themes.R")
require(pomp)
require(panelPomp)
select <- dplyr::select
summarize <- dplyr::summarize
test_data_filename <- "H1_data_serology_swab.rda"
load(test_data_filename)
serology_sub <- serology_sub_H1
serology_sub$age = NA
for( i in c(1:nrow(serology_sub))){
  serology_sub[i,]$age = demog[demog$memberID == serology_sub[i,]$memberID,]$age_at_recruitment
}

# Calculate the time between first titer and swab, and between swab and second titer
sub <- serology_sub%>% 
  mutate(Child  = as.factor(age < 15)) %>% 
  group_by(memberID,swab_date,Child,age,min_date,max_date) %>% 
  summarize(delta = log_titer_trans(max(value)) - log_titer_trans(min(value)),
            pre = log_titer_trans(min(value)),
            dt1 = as.numeric(max(swab_date) - min(date)),
            dt2 = as.numeric(max(date) - max(swab_date)),
            dt_total = as.numeric(max(date) - min(date))) 

p1 <- ggplot(sub %>% filter(dt_total < 365), aes(x = pre, y = delta)) + 
  geom_point() + 
  stat_smooth(method = "lm", color = "red", se = F) + 
  geom_jitter(height = .1, width = .1) + 
  xlab("Log pre-infection titer") + 
  ylab(expression(Delta~titer)) + 
  plot_themes

p2 <- ggplot(sub %>% filter(dt_total < 365), aes(x = pre, y = delta)) + 
  geom_point() + 
  facet_wrap(~Child) + 
  stat_smooth(method = "lm", color = "red", se = F) + 
  geom_jitter(height = .1, width = .1) + 
  xlab("Log pre-infection titer") + 
  ylab(expression(Delta~titer)) + 
  plot_themes

p_composite <- plot_grid( p1, p2, ncol = 1)
save_plot("delta_vs_pre.pdf", p_composite, base_width = 6, base_height = 4)
