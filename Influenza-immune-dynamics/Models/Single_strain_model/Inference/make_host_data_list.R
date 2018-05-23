require(ggplot2)
require(MASS)
require(RSQLite)
require(reshape2)
require(dplyr)
source("../Utility_scripts/model_functions.R")
source("../../../Imprinting/imprinting_functions.R")
require(panelPomp)
select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains 

## Set up input and output files  ## ----------------------------------------------------------------------------------------------------------------------------------
dbFilename <- "../../../Data/Data.sqlite"
test_subtype = "pH1N1"
host_data_filename = paste0("pomp_data_",test_subtype, ".rda")
calculate_imprinting = T

#Read in data  
db <- dbConnect(SQLite(), dbFilename)
serology <- dbReadTable(db, "serology_full_model") %>% 
  filter(subtype == test_subtype) %>% 
  mutate(date = as.Date(date, origin = "1970-1-1"))
demography <- dbReadTable(db, "demography_full_model") %>% 
  mutate(birthdate = as.Date(birthdate, origin = "1970-1-1"),
         visitdate = as.Date(visitdate, origin = "1970-1-1"))
dbDisconnect(db)

## Imprinting data ## --------------------------------------------------------------------------------------------------
if(calculate_imprinting){
  duration_maternal_Abs = 6 # Months
  max_age_first_exposure = 12
  MONTHS_PER_YEAR = 12
  smooth_frequency_data = F
  
  pandemic_data <- read.csv("../../../Imprinting/frequency_data/Historic_Flu_Pandemics.csv") %>% rename(year = FirstYearOfSeason,
                                                                                                        frac_h3n2 = H3N2_fraction,
                                                                                                        frac_h1n1 = H1N1_fraction,
                                                                                                        frac_h2n2 = H2N2_fraction) %>%    select(-Source) %>% 
    filter(year < 1968) %>% 
    arrange(year) %>% 
    select(year,frac_h1n1,frac_h3n2, frac_h2n2)
  pandemic_fractions <- pandemic_data
  
  GISAID_data = read.csv("../../../Imprinting/frequency_data/GISAID_fractions.csv") %>% 
    select(-X) %>% 
    filter(year < 1997)
  
  flunet_data = read.csv("../../../Imprinting/frequency_data/FluNetInteractiveReport.csv") %>% mutate(frac_h1n1 = (AH1/INF_A))
  flunet_yearly <- flunet_data %>% group_by(Year) %>% summarize(h1 = sum(AH1, AH1N12009, na.rm=T), h3 = sum(AH3,na.rm=T), all_A = sum(INF_A,na.rm=T)) %>% 
    mutate(frac_h1n1 = h1/all_A, frac_h3n2 = h3/all_A, frac_h2n2 = 0) %>% 
    rename(year = Year) %>% 
    select(year, frac_h1n1, frac_h3n2,frac_h2n2) %>% 
    filter(year %in% c(1997:2012))
  
  relative_incidence <- rbind(pandemic_fractions, GISAID_data,flunet_yearly)
  imprinting_data <- lapply(unique(serology$memberID), prob_primary_exposure, demog_data = demography, sero_data = serology)
  df_imprinting <- as.data.frame(do.call("rbind",imprinting_data)) %>% mutate(memberID = unique(serology$memberID))
  save(df_imprinting, file = imprinting_data_filename)
}
if(!calculate_imprinting){
  load(imprinting_data_filename)
}

## Make list of host data -----------------------------------------------------------------------------------------------------
strains = subtype
host_data_list <- lapply(test_ids, 
                         FUN = make_pomp_data_for_ind, 
                         demographic_data = demog,
                         strain_vec = strains,
                         imprinting_data = df_imprinting,
                         simulate_times = F, # Do not simulate additional visits for individuals
                         n_vis_sim = NA) # Number of observations per individual IF simulating data

save(host_data_list, test_ids, file = host_data_filename)