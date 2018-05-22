## Extract and format data 
library(RSQLite)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lubridate)
source("../Utility_scripts/data_formatting_functions.R")
source("imprinting_functions.R")
textSize = 12
source("../Utility_scripts/plot_themes.R")
select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise

## Get the participant data ## ----------------------------------------------------------------------------------------
dbFilename <- "../Data/hongkongflu.sqlite"
#source("../Raw_data_analysis/extract_relevant_data.R")
load("all_data.rda")
## ----------------------------------------------------------------------------------------

duration_maternal_Abs = 6 # in months 
max_age_first_exposure = 12
MONTHS_PER_YEAR = 12
smooth_frequency_data = F
make_yearly_dfs = F
test_years <- c(2009:2014)
df_filenames <- paste0("df_imprinting_",test_years, ".rda")

## Smooth the monthly frequency data ## ----------------------------------------------------------------------------------
if(smooth_frequency_data){
  strains <- c("h1n1","h3n2","h2n2")
  for( i in c(1:length(strains))){
    this_strain <- strains[i]
    dat <- relative_incidence %>%  select(date, contains(this_strain))
    names(dat) <- c("x","y")
    dat$x <- as.numeric(dat$x)
    smoo <- with(dat[!is.na(dat$y),],smooth.spline(x,y))
    result <- with(dat,predict(smoo,x[is.na(y)]))
    dat[is.na(dat$y),] <- result
    relative_incidence[,grep(this_strain,names(relative_incidence))] <- dat$y
  }
}
if(make_yearly_dfs){
  pandemic_data <- read.csv("../Imprinting/frequency_data/Historic_Flu_Pandemics.csv") %>% rename(year = FirstYearOfSeason,
                                                                                                  frac_h3n2 = H3N2_fraction,
                                                                                                  frac_h1n1 = H1N1_fraction,
                                                                                                  frac_h2n2 = H2N2_fraction) %>% 
    select(-Source) %>% 
    filter(year < 1968) %>% 
    arrange(year) %>% 
    select(year,frac_h1n1,frac_h3n2, frac_h2n2)
  pandemic_fractions <- pandemic_data
  
  GISAID_fractions = read.csv("./frequency_data/GISAID_fractions.csv") %>% 
    select(year, frac_h1n1, frac_h3n2,frac_h2n2) 
  
  flunet_data = read.csv("../Imprinting/frequency_data/FluNetInteractiveReport.csv") 
  flunet_yearly <- flunet_data %>% group_by(Year) %>% summarize(h1 = sum(AH1, AH1N12009, na.rm=T), h3 = sum(AH3,na.rm=T)) %>% #, all_A = sum(INF_A,na.rm=T)) %>% 
    mutate(frac_h1n1 = h1/(h1 + h3), frac_h3n2 = h3/(h1+h3), frac_h2n2 = 0) %>% 
    rename(year = Year) %>% 
    select(year, frac_h1n1, frac_h3n2,frac_h2n2) %>% 
    filter(year %in% c(1997:2012))
  
  relative_incidence <- rbind(pandemic_fractions, genbank_fractions,flunet_yearly)

  for( i in c(1:length(test_years))){
    this_year = test_years[i]
    cat("this year is ",this_year, "\n")
    month = 1
    imprinting_data <- lapply(test_ids, 
                              prob_primary_exposure, 
                              demog_data = demog, 
                              sero_data = serology,
                              test_year = this_year,
                              test_month = 1
                             )
    df_imprinting <- as.data.frame(do.call("rbind",imprinting_data)) %>% mutate(memberID = test_ids)
    filename <- paste0("df_imprinting_", this_year, ".rda")
    save(df_imprinting, file = filename)
  }
}

df_all_years <- data.frame()
for(i in c(1:length(test_years))){
  test_year = test_years[i]
  load(df_filenames[i])
  df_imp <- df_imprinting
  df_imp$test_year = test_year
  df_imp$DOB <- NA
  for( i in c(1:nrow(df_imprinting))){
    df_imp[i,]$DOB <- demog %>% filter(memberID == df_imp[i,]$memberID) %>% select(birthdate) %>% unlist
  }
  df_imp$DOB <- as.Date(df_imp$DOB, origin = "1970-1-1")
  df_imp <- df_imp %>% mutate(age = test_year - year(df_imp$DOB))
  df_all_years <- rbind(df_all_years, df_imp)
}

df_all_years_sum <- df_all_years %>% group_by(test_year, age) %>% summarize(H1N1 = mean(h1n1),
                                                                            H3N2 = mean(h3n2),
                                                                            H2N2 = mean(h2n2),
                                                                            Naive = mean(p_naive)
                                                                            ) %>% as.data.frame()
df_age_delta <- data.frame()
ages <- unique(df_all_years_sum$age)

for(i in c(1:length(ages))){
  this_age = ages[i]
  df_sub <- df_all_years_sum %>% filter(age == this_age & !is.na(H1N1) & !is.na(H3N2) & !is.na(H2N2) & !is.na(Naive))
  df_max<- unlist(df_sub[df_sub$test_year == max(df_sub$test_year),])
  df_min <- unlist(df_sub[df_sub$test_year == min(df_sub$test_year),])
  delta = df_max[-c(1:2)] - df_min[-c(1:2)]
  df_age <- data.frame(c(age = this_age,
                       as.data.frame(t(delta))))
  df_age_delta <- rbind(df_age_delta, df_age)
}

age_pop_df <- df_all_years %>% select(memberID, test_year, age) %>% group_by(test_year, age) %>% summarize(n = n_distinct(memberID))
age_pop_sum <- age_pop_df %>% group_by(age) %>% summarize(mean_n = mean(n))
dfm_delta <- melt(df_age_delta, id.vars = "age") %>% 
  mutate(value = round(value,2))


p_delta <- ggplot(dfm_delta, aes( x= age, y = value)) +
  geom_line(aes(color = variable)) +
  xlab("Age") +
  xlim(6,70) +
  ylab(expression(Delta*" probability of primary infection")) +
  labs("") +
  plot_themes

p_delta <- p_delta + geom_line(data = age_pop_sum, aes(x = age, y = mean_n/40),linetype =2) +
  scale_y_continuous(sec.axis = sec_axis(~.*40, name = "Mean number of individuals"))

snapshot_sub <- df_all_years_sum %>% filter(test_year == 2009) %>% 
  select(H1N1,H3N2,H2N2,Naive,age)
dfm <- melt(snapshot_sub, id.vars = "age")%>% 
  mutate(value = round(value,2))

p_snapshot <- ggplot(dfm, aes(x = as.factor(age), y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity") +
  ylim(0,1) + 
  plot_themes + 
  xlab("")


if(save_plots){
  save_plot("change_in_imprinting_probs.pdf", p_delta, base_width = 10, base_height = 8)
  save_plot("snapshot_imprinting.pdf", p_snapshot, base_width = 10, base_height = 8)
}
