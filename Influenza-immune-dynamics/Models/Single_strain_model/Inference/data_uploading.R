
load("../../Short_term_sub_model/Raw_data_analysis/H1_data_serology_swab.rda")
dbFilename = "../../../Data/Data.sqlite"
db <- dbConnect(SQLite(), dbFilename)
#dbWriteTable(db, "serology_full_model", serology , overwrite = T)
#dbWriteTable(db, "demography_full_model", demog , overwrite = T)
#dbWriteTable(db, "data_sub_model_H3", serology_sub_H3, overwrite = T)
dbDisconnect(db)


demography_full <- dbReadTable(db, "demography")
flu_pop_full <- dbReadTable(db, "flu_population")
participants_full <- dbReadTable(db, "participantlist")
serology_full <- dbReadTable(db, "serology_full_model") %>% 
  mutate(date = as.Date(date, origin = "1970-1-1"))
serology_visit_dates_full <- dbReadTable(db,"serology_visit_dates")
swabs_full <- dbReadTable(db,"swabs")
dbDisconnect(db)
