## Extract and format data 
library(RSQLite)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
source("../Utility_scripts/data_formatting_functions.R")
source("../Imprinting/imprinting_functions.R")
textSize = 12
source("../Utility_scripts/plot_themes.R")
select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise

## Get the participant data ## ----------------------------------------------------------------------------------------
dbFilename <- "../Data/hongkongflu.sqlite"
include_pilot = F
pilot_data_filename <- "../Raw_data_analysis/extract_pilot_data_pH1.R"
source("../Raw_data_analysis/extract_relevant_data.R")

## Imprinting data ## --------------------------------------------------------------------------------------------------
calculate_imprinting = F

if(calculate_imprinting){
  duration_maternal_Abs = 6 # in months 
  max_age_first_exposure = 12
  MONTHS_PER_YEAR = 12
  smooth_frequency_data = F
  
  pandemic_data <- read.csv("../Imprinting/frequency_data/Historic_Flu_Pandemics.csv") %>% rename(year = FirstYearOfSeason,
                                                                                      frac_h3n2 = H3N2_fraction,
                                                                                      frac_h1n1 = H1N1_fraction,
                                                                                      frac_h2n2 = H2N2_fraction) %>% 
    select(-Source) %>% 
    filter(year < 1968) %>% 
    arrange(year) %>% 
    select(year,frac_h1n1,frac_h3n2, frac_h2n2)
  pandemic_fractions <- pandemic_data
  
  genbank_data_China = read.csv("../Imprinting/frequency_data/subtype_fraction_6812_HK_China_byYear.csv") 
  genbank_data_global = read.csv("../Imprinting/frequency_data/subtype_fraction_6812_Global_byYear.csv")
  genbank_fractions <- genbank_data_global %>% mutate(a = sum(h3n2,h1n1, na.rm=T)) %>% 
    mutate(frac_h1n1 = h1n1/a, frac_h3n2 = h3n2/a) %>%
    select(year,frac_h1n1,frac_h3n2) %>% 
    mutate(frac_h2n2 = 0) %>% 
    filter(year < 1997)
  
  flunet_data = read.csv("../Imprinting/frequency_data/FluNetInteractiveReport.csv") %>% mutate(frac_h1n1 = (AH1/INF_A))
  flunet_yearly <- flunet_data %>% group_by(Year) %>% summarize(h1 = sum(AH1, AH1N12009, na.rm=T), h3 = sum(AH3,na.rm=T), all_A = sum(INF_A,na.rm=T)) %>% 
    mutate(frac_h1n1 = h1/all_A, frac_h3n2 = h3/all_A, frac_h2n2 = 0) %>% 
    rename(year = Year) %>% 
    select(year, frac_h1n1, frac_h3n2,frac_h2n2) %>% 
    filter(year %in% c(1997:2012))
  
  relative_incidence <- rbind(pandemic_fractions, genbank_fractions,flunet_yearly) #%>% select(year, month, contains("frac"))
  #relative_incidence$date = as.Date( paste( relative_incidence$year , relative_incidence$month ,"1", sep = "-" )  , format = "%Y-%m-%d")
  
  ## Smooth the monthly frequency data ## ----------------------------------------------------------------------------------
  if(smooth_frequency_data){
    strains <- c("h1n1","h3n2","b","h2n2")
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
  
  imprinting_data <- lapply(test_ids, prob_primary_exposure, demog_data = demog, sero_data = serology)
  df_imprinting <- as.data.frame(do.call("rbind",imprinting_data)) %>% mutate(memberID = test_ids)
  save(df_imprinting, file = "../Imprinting/imprinting_data_subjects_full_data.rda")
}
if(!calculate_imprinting){
  load("../Imprinting/imprinting_data_subjects_full_data.rda")
}
##----------------------------------------------------------------------------------------------------------------------

## Generate a list of host data objects and save output
strains = unique(serology$name)[c(3,5)]
host_data_list <- lapply(test_ids, 
                         FUN = make_pomp_data_for_ind, 
                         demographic_data = demog,
                         strain_vec = strains,
                         imprinting_data = df_imprinting,
                         simulate_times = T, 
                         n_vis_sim = 150)
save(host_data_list, test_ids, file = "pomp_data_H1_H3_kids_sim.rda")
##----------------------------------------------------------------------------------------------------------------------
