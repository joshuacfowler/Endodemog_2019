## Authors: Josh and Tom	## Grass endophyte population model
## Purpose: Create a script that imports Endodemog data, perform all raw data manipulation to set up flowering tiller and seed production info,	
## and create an .RData object that can be loaded for analysis	
## Last Update: 11/28/2018
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)

# This is the main updated data sheet that we will merge with the flowering and seed data
LTREB_endodemog <- 
  read.csv(file = "endo_demog_long.csv")


str(LTREB_endodemog)
dim(LTREB_endodemog)


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(species_0 = (recode_factor(species,                   
                                    "AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                                    "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                                    "POSY" = 7))) %>%        
  mutate(species_index = as.integer(species_0)) %>% 
  mutate(year_t_index = as.integer(recode_factor(year_t, 
                                                 '2007' = 1, '2008' = 2, '2009' = 3, 
                                                 '2010' = 4, '2011' = 5, '2012' = 6, 
                                                 '2013' = 7, '2014' = 8, '2015' = 9, 
                                                 '2016' = 10, '2017' = 11))) %>%             
  mutate(year_t1_index = as.integer(recode_factor(year_t1, 
                                                  '2008' = 2, '2009' = 3, '2010' = 4, 
                                                  '2011' = 5, '2012' = 6, '2013' = 7, 
                                                  '2014' = 8, '2015' = 9, '2016' = 10, 
                                                  '2017' = 11, '2018' = 12))) %>%               
  mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, 
                                          origin == "R" ~ 1, 
                                          origin != "R" | origin != "O" ~ 1))) %>%   
  mutate(plot_fixed = (case_when(plot != "R" ~ plot, 
                                 plot == "R" ~ origin))) %>%                         
  mutate(plot_index = as.integer(as.factor(plot_fixed))) %>%                         
  mutate(endo_01 = case_when(endo == "0" | endo == "minus" ~ 0,
                             endo == "1"| endo =="plus" ~ 1)) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))                              

dim(LTREB_data)
unique(LTREB_data$surv_t1)


# read in data from POAL
POAL_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALnew complete with 2016 data.xlsx", sheet = "POAL")
POAL_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALnew complete with 2016 data.xlsx", sheet = "POAL (NEW) recruits")
POAL_data_old <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALold complete with 2016 data.xlsx", sheet = "POAL (OLD)")
POAL_data_old_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALold complete with 2016 data.xlsx", sheet = "POAL (OLD) recruits")

# read in data from POSY
POSY_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYnew complete with 2016 data.xlsx", sheet = "POSY")
POSY_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYnew complete with 2016 data.xlsx", sheet = "POSY (NEW) recruits")
POSY_data_old <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYold complete with 2016 data.xlsx", sheet = "POSY(Old)")
POSY_data_old_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYold complete with 2016 data.xlsx", sheet = "POSY(Old)recruits")

# read in data from LOAR
LOAR_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_to2016_complete.xlsx", sheet = "LOAR")
LOAR_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_to2016_complete.xlsx", sheet = "LOAR recruits")
LOAR_data_r <- LOAR_data_r%>% # assign a Tag name for recruits and endophyte status based on plot and id
  mutate(tag = paste(Plot, `Recruit ID`, sep = "-")) %>% 
  mutate(Endo = LOAR_data$Endo[match(LOAR_data_r$Plot, LOAR_data$PLOT)])

# Read in data from FESU
FESU_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU")
FESU_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU recruits")

# Read in data from ELVI
ELVI_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI(IN) originals up to 2016.xlsx", sheet = "ELVI")
ELVI_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVI recruits")
ELVI_data_r <- ELVI_data_r %>% 
  mutate(tag = paste(Plot, RecruitNo, sep = "-")) 
  

# Read in data from ELRI
ELRI_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI")
ELRI_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI recruits")
ELRI_data_r <- ELRI_data_r %>% 
  mutate(tag = paste(PLOT, RecruitNo, sep = "-")) 
# Read in data from POSY
AGPE_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/AGPE2016_final.xlsx", sheet = "AGPE")
AGPE_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/AGPE2016_final.xlsx", sheet = "AGPE recruits")

# A list of the columns that are redundant for now, or will be used in our seed estimates
 # pmain <- POAL_data %>%  
 #  select(-survive4may08,-`Surv4/11`, -Aphids1, -aphids2, -aphid3, -contains("SB"),-notes2,
 #         -notes2__1, -notes3, -notes4, -notes5,-`notes5/4/2008`,-notes6, -notes7, 
 #         -notes8, -notes9, -othernotes3, -data, -`Planting Notes`, -TAG, -Endocheck, 
 #         -EndoDateCheck, -EndoDateCheck_Day, -EndoDateCheck_Month, -EndoDateCheck_Year,
 #         -`Coll Date`, -TotTillers11sep08, -tilleradjust, -VisualEST1, -endoyr, 
 #         -seed2surv, -seed3surv, -seed4surv, -actseed2, -contains("Est"), -contains("Hbv"), 
 #         -contains("Lvs"), -contains("Infl"), -contains("CB"), -contains("Prop"))

# Combining measurements across years for the “New” POAL data ------------------

## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

psurv <- POAL_data %>%
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("survive1", "survive2", "Survive3", "Survive4", 
                       "Survive5", "survive6", "survive7", "survive8", 
                       "survive9"),
       value.name = "surv") 
psurv$year<- ifelse(psurv$variable == "survive1", 2008, ifelse(psurv$variable  == "survive2", 2009, ifelse(psurv$variable  == "Survive3", 2010, ifelse(psurv$variable  == "Survive4", 2011, ifelse(psurv$variable  == "Survive5", 2012, ifelse(psurv$variable  == "survive6", 2013,ifelse(psurv$variable == "survive7", 2014,ifelse(psurv$variable == "survive8", 2015,ifelse(psurv$variable  == "survive9", 2016, NA)))))))))
# View(psurv)

pgrow <- POAL_data %>% 
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("Tottillers1", "TotTillers2", "TotTillers3",
                       "TotTillers4", "TotTillers5", "TotTillers6", 
                       "TotTillers7", "TotTillers8", "TotTillers9"), 
       value.name = "size") 
pgrow$year<- ifelse(pgrow$variable == "Tottillers1", 2008, ifelse(pgrow$variable  == "TotTillers2", 2009, ifelse(pgrow$variable  == "TotTillers3", 2010, ifelse(pgrow$variable  == "TotTillers4", 2011, ifelse(pgrow$variable  == "TotTillers5", 2012, ifelse(pgrow$variable  == "TotTillers6", 2013, ifelse(pgrow$variable == "TotTillers7", 2014, ifelse(pgrow$variable == "TotTillers8", 2015, ifelse(pgrow$variable  == "TotTillers9", 2016, NA)))))))))
# View(pgrow)

pflw <- POAL_data %>% 
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
pflw$year<- ifelse(pflw$variable == "Flwtillers1", 2008, ifelse(pflw$variable  == "FlwTillers2", 2009, ifelse(pflw$variable  == "FlwTillers3", 2010, ifelse(pflw$variable  == "FlwTillers4", 2011, ifelse(pflw$variable  == "FlwTillers5", 2012, ifelse(pflw$variable  == "FlwTillers6", 2013,ifelse(pflw$variable == "FlwTillers7", 2014,ifelse(pflw$variable == "FlwTillers8", 2015,ifelse(pflw$variable  == "FlwTillers9", 2016, NA)))))))))
# View(pflw)

pmerge_sg <- merge(psurv, pgrow, by = c( "plot", "pos", "tag", "Endo", 
                                              "Loc'n", "Birth Year", "TRT",
                                              "Plant", "year"))
# View(pmerge_sg)

pmerge_sgf <- merge(pmerge_sg, pflw, by = c( "plot", "pos", "tag", "Endo", 
                                             "Loc'n", "Birth Year", "TRT",
                                             "Plant", "year"))
# View(pmerge_sgf)

# getting a dataframe with t and t_1
pmerge_t1 <-pmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(pmerge_t1)

pmerge_t <-pmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(pmerge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
pmerge <- pmerge_t1 %>% 
  full_join(pmerge_t, by = c("plot", "pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(pmerge)


# Combining measurements across years for the "New" POAL recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
rsurv <- POAL_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive10 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
rsurv$year<- ifelse(rsurv$variable == "survive10", 2010, ifelse(rsurv$variable == "Survive11", 2011, ifelse(rsurv$variable  == "Survive12", 2012, ifelse(rsurv$variable  == "Survive13", 2013, ifelse(rsurv$variable  == "Survive14", 2014, ifelse(rsurv$variable  == "Survive15", 2015, ifelse(rsurv$variable  == "Survive16", 2016, NA)))))))
# View(rsurv)

rgrow <- POAL_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
rgrow$year<- ifelse(rgrow$variable == "TOTtiller10", 2010, ifelse(rgrow$variable == "TOTtiller11", 2011, ifelse(rgrow$variable  == "TOTtiller12", 2012, ifelse(rgrow$variable  == "TOTtiller13", 2013, ifelse(rgrow$variable  == "TOTtiller14", 2014, ifelse(rgrow$variable  == "TOTtiller15", 2015, ifelse(rgrow$variable  == "TOTtiller16", 2016, NA)))))))
# View(rgrow)

rflw <- POAL_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
rflw$year<- ifelse(rflw$variable == "FLWtiller10", 2010, ifelse(rflw$variable == "FLWtiller11", 2011, ifelse(rflw$variable  == "FLWtiller12", 2012, ifelse(rflw$variable  == "FLWtiller13", 2013, ifelse(rflw$variable  == "FLWtiller14", 2014, ifelse(rflw$variable  == "FLWtiller15", 2015, ifelse(rflw$variable  == "FLWtiller16", 2016, NA)))))))
# View(rflw)

rmerge_sg <- merge(rsurv, rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(rmerge_sg)

rmerge_sgf <- merge(rmerge_sg, rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(rmerge_sgf)

## getting a dataframe with time t and t_1
rmerge_t1 <-rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(rmerge_t1)

rmerge_t <-rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(rmerge_t)

rmerge <- rmerge_t1 %>% 
  full_join(rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(rmerge)




# Combining measurements across years for the “Old” POAL data ------------------


## Combining data across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
poldsurv <- POAL_data_old %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"),
       measure.var = c("survive1", "Survive3", "survive4", 
                       "survive5", "survive6", "survive7", 
                       "survive8", "survive9", "survive10"),
       value.name = "surv") 
poldsurv$year<- ifelse(poldsurv$variable == "survive1", 2008, ifelse(poldsurv$variable  == "Survive3", 2009, ifelse(poldsurv$variable  == "survive4", 2010, ifelse(poldsurv$variable  == "survive5", 2011, ifelse(poldsurv$variable  == "Survive6", 2012, ifelse(poldsurv$variable  == "survive7", 2013,ifelse(poldsurv$variable == "survive8", 2014,ifelse(poldsurv$variable == "survive9", 2015,ifelse(poldsurv$variable  == "survive10", 2016, NA)))))))))
# View(poldsurv)

poldgrow <- POAL_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("TotTillers1", "TotTillers3","TotTillers4", 
                       "TotTillers5", "TotTillers6","TotTillers7", 
                       "TotTillers8", "TotTillers9", "TotTillers10"), 
       value.name = "size") 
poldgrow$year<- ifelse(poldgrow$variable == "TotTillers1", 2008, ifelse(poldgrow$variable  == "TotTillers3", 2009, ifelse(poldgrow$variable  == "TotTillers4", 2010, ifelse(poldgrow$variable  == "TotTillers5", 2011, ifelse(poldgrow$variable  == "TotTillers6", 2012, ifelse(poldgrow$variable  == "TotTillers7", 2013, ifelse(poldgrow$variable == "TotTillers8", 2014, ifelse(poldgrow$variable == "TotTillers9", 2015, ifelse(poldgrow$variable  == "TotTillers10", 2016, NA)))))))))
 # View(poldgrow)

poldflw <- POAL_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n",
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("Flwtillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
poldflw$year<- ifelse(poldflw$variable == "Flwtillers1", 2008, ifelse(poldflw$variable  == "FlwTillers3", 2009, ifelse(poldflw$variable  == "FlwTillers4", 2010, ifelse(poldflw$variable  == "FlwTillers5", 2011, ifelse(poldflw$variable  == "FlwTillers6", 2012, ifelse(poldflw$variable  == "FlwTillers7", 2013,ifelse(poldflw$variable == "FlwTillers8", 2014,ifelse(poldflw$variable == "FlwTillers9", 2015,ifelse(poldflw$variable  == "FlwTillers10", 2016, NA)))))))))
# View(poldflw)


poldmerge_sg <- merge(poldsurv, poldgrow, by = c( "plot","pos", "tag", "Endo", 
                                         "Loc'n", "Birth Year", "TRT",
                                         "Plant", "year"))
# View(poldmerge_sg)

poldmerge_sgf <- merge(poldmerge_sg, poldflw, by = c( "plot","pos", "tag", "Endo", 
                                             "Loc'n", "Birth Year", "TRT",
                                             "Plant", "year"))
# View(poldmerge_sgf)

# getting a dataframe with t and t_1
poldmerge_t1 <-poldmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(pmerge_t1)

poldmerge_t <-poldmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(pmerge_t)

poldmerge <- poldmerge_t1 %>% 
  full_join(poldmerge_t, by = c("plot","pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
 # View(poldmerge)


# Combining measurements across years for the “Old” POAL recruits data --------

## Combining data for recruits across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
roldsurv <- POAL_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive09 = "NA") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("survive09", "Survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
roldsurv$year<- ifelse(roldsurv$variable == "survive09", 2009, ifelse(roldsurv$variable == "Survive10", 2010, ifelse(roldsurv$variable == "Survive11", 2011, ifelse(roldsurv$variable  == "Survive12", 2012, ifelse(roldsurv$variable  == "Survive13", 2013, ifelse(roldsurv$variable  == "Survive14", 2014, ifelse(roldsurv$variable  == "Survive15", 2015, ifelse(roldsurv$variable  == "Survive16", 2016, NA))))))))
# View(roldsurv)


roldgrow <- POAL_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller09", "TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
roldgrow$year<- ifelse(roldgrow$variable == "TOTtiller09", 2009, ifelse(roldgrow$variable == "TOTtiller10", 2010, ifelse(roldgrow$variable == "TOTtiller11", 2011, ifelse(roldgrow$variable  == "TOTtiller12", 2012, ifelse(roldgrow$variable  == "TOTtiller13", 2013, ifelse(roldgrow$variable  == "TOTtiller14", 2014, ifelse(roldgrow$variable  == "TOTtiller15", 2015, ifelse(roldgrow$variable  == "TOTtiller16", 2016, NA))))))))
# View(roldgrow)

roldflw <- POAL_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWTiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
roldflw$year<- ifelse(roldflw$variable == "FLWTiller09", 2009, ifelse(roldflw$variable == "FLWtiller10", 2010, ifelse(roldflw$variable == "FLWtiller11", 2011, ifelse(roldflw$variable  == "FLWtiller12", 2012, ifelse(roldflw$variable  == "FLWtiller13", 2013, ifelse(roldflw$variable  == "FLWtiller14", 2014, ifelse(roldflw$variable  == "FLWtiller15", 2015, ifelse(roldflw$variable  == "FLWtiller16", 2016, NA))))))))
# View(roldflw)

roldmerge_sg <- merge(roldsurv, roldgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(roldmerge_sg)

roldmerge_sgf <- merge(roldmerge_sg, roldflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(roldmerge_sgf)

## getting a dataframe with time t and t_1
roldmerge_t1 <-roldmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(roldmerge_t1)

roldmerge_t <-roldmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(roldmerge_t)

roldmerge <- roldmerge_t1 %>% 
  full_join(roldmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(roldmerge)










# Combining the old and new and original and recruit POAL dataframes ---------
pmerge <- pmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                     "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                     "year_t", "size_t", "flw_t")]
poldmerge <- poldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "size_t", "flw_t")]
rmerge <- rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                   "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                   "year_t", "size_t", "flw_t")]
roldmerge <- roldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "size_t", "flw_t")]
# View(pmerge)
# View(poldmerge)
# View(rmerge)
# View(roldmerge)



POAL <- pmerge %>% 
  rbind(poldmerge) %>% 
  rbind(roldmerge) %>% 
  rbind(rmerge) %>% 
  mutate(species = "POAL") %>% 
  mutate(surv_t1 = ifelse(surv_t1 == 'NA', NA, surv_t1))
POAL <- POAL[!(is.na(POAL$surv_t1)),]

View(POAL)









# Combining measurements across years for the “New” POSY data -------------

## recoding for the year of measurement
## merging these measurements into one dataframe

po_surv <- POSY_data %>%
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("survive1", "survive2", "Survive3", "survive4", 
                       "survive5", "survive6", "survive7", "survive8"),
       value.name = "surv") 
po_surv$year<- ifelse(po_surv$variable == "survive1", 2008, ifelse(po_surv$variable  == "survive2", 2009, ifelse(po_surv$variable  == "Survive3", 2010, ifelse(po_surv$variable  == "survive4", 2011, ifelse(po_surv$variable  == "survive5", 2012, ifelse(po_surv$variable  == "survive6", 2013,ifelse(po_surv$variable == "survive7", 2014,ifelse(po_surv$variable == "survive8", 2015,ifelse(po_surv$variable  == "survive9", 2016, NA)))))))))
# View(po_surv)

po_grow <- POSY_data %>% 
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("Tottillers1", "TotTillers2", "TotTillers3",
                       "TotTillers4", "TotTillers5", "TotTillers6", 
                       "TotTillers7", "TotTillers8"), 
       value.name = "size") 
po_grow$year<- ifelse(po_grow$variable == "Tottillers1", 2008, ifelse(po_grow$variable  == "TotTillers2", 2009, ifelse(po_grow$variable  == "TotTillers3", 2010, ifelse(po_grow$variable  == "TotTillers4", 2011, ifelse(po_grow$variable  == "TotTillers5", 2012, ifelse(po_grow$variable  == "TotTillers6", 2013, ifelse(po_grow$variable == "TotTillers7", 2014, ifelse(po_grow$variable == "TotTillers8", 2015, ifelse(po_grow$variable  == "TotTillers9", 2016, NA)))))))))
# View(po_grow)

po_flw <- POSY_data %>% 
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8"), 
       value.name = "flw") 
po_flw$year<- ifelse(po_flw$variable == "Flwtillers1", 2008, ifelse(po_flw$variable  == "FlwTillers2", 2009, ifelse(po_flw$variable  == "FlwTillers3", 2010, ifelse(po_flw$variable  == "FlwTillers4", 2011, ifelse(po_flw$variable  == "FlwTillers5", 2012, ifelse(po_flw$variable  == "FlwTillers6", 2013,ifelse(po_flw$variable == "FlwTillers7", 2014,ifelse(po_flw$variable == "FlwTillers8", 2015,ifelse(po_flw$variable  == "FlwTillers9", 2016, NA)))))))))
# View(po_flw)

po_merge_sg <- merge(po_surv, po_grow, by = c( "plot", "pos", "tag", "Endo", 
                                         "Loc'n", "Birth Year", "TRT",
                                         "Plant", "year"))
# View(po_merge_sg)

po_merge_sgf <- merge(po_merge_sg, po_flw, by = c( "plot", "pos", "tag", "Endo", 
                                             "Loc'n", "Birth Year", "TRT",
                                             "Plant", "year"))
# View(po_merge_sgf)

# getting a dataframe with t and t_1
po_merge_t1 <-po_merge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(po_merge_t1)

po_merge_t <-po_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(po_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits

po_merge <- po_merge_t1 %>% 
  full_join(po_merge_t, by = c("plot", "pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(po_merge)


# Combining measurements across years for the "New" POSY recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
po_rsurv <- POSY_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive10 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
po_rsurv$year<- ifelse(po_rsurv$variable == "survive10", 2010, ifelse(po_rsurv$variable == "Survive11", 2011, ifelse(po_rsurv$variable  == "Survive12", 2012, ifelse(po_rsurv$variable  == "Survive13", 2013, ifelse(po_rsurv$variable  == "Survive14", 2014, ifelse(po_rsurv$variable  == "Survive15", 2015, ifelse(po_rsurv$variable  == "Survive16", 2016, NA)))))))
# View(po_rsurv)

po_rgrow <- POSY_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
po_rgrow$year<- ifelse(po_rgrow$variable == "TOTtiller10", 2010, ifelse(po_rgrow$variable == "TOTtiller11", 2011, ifelse(po_rgrow$variable  == "TOTtiller12", 2012, ifelse(po_rgrow$variable  == "TOTtiller13", 2013, ifelse(po_rgrow$variable  == "TOTtiller14", 2014, ifelse(po_rgrow$variable  == "TOTtiller15", 2015, ifelse(po_rgrow$variable  == "TOTtiller16", 2016, NA)))))))
# View(po_rgrow)

po_rflw <- POSY_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_rflw$year<- ifelse(po_rflw$variable == "FLWtiller10", 2010, ifelse(po_rflw$variable == "FLWtiller11", 2011, ifelse(po_rflw$variable  == "FLWtiller12", 2012, ifelse(po_rflw$variable  == "FLWtiller13", 2013, ifelse(po_rflw$variable  == "FLWtiller14", 2014, ifelse(po_rflw$variable  == "FLWtiller15", 2015, ifelse(po_rflw$variable  == "FLWtiller16", 2016, NA)))))))
# View(po_rflw)

po_rmerge_sg <- merge(po_rsurv, po_rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_rmerge_sg)

po_rmerge_sgf <- merge(po_rmerge_sg, po_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_rmerge_sgf)

## getting a dataframe with time t and t_1
po_rmerge_t1 <-po_rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(po_rmerge_t1)

po_rmerge_t <-po_rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(po_rmerge_t)

po_rmerge <- po_rmerge_t1 %>% 
  full_join(po_rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(po_rmerge)




# Combining measurements across years for the “Old” POSY data ------------------


## Combining data across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
po_oldsurv <- POSY_data_old %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"),
       measure.var = c("survive1", "survive3", "survive4", 
                       "survive5", "Survive6", "survive7", 
                       "survive8", "survive9", "survive10"),
       value.name = "surv") 
po_oldsurv$year<- ifelse(po_oldsurv$variable == "survive1", 2008, ifelse(po_oldsurv$variable  == "survive3", 2009, ifelse(po_oldsurv$variable  == "survive4", 2010, ifelse(po_oldsurv$variable  == "survive5", 2011, ifelse(po_oldsurv$variable  == "Survive6", 2012, ifelse(po_oldsurv$variable  == "survive7", 2013,ifelse(po_oldsurv$variable == "survive8", 2014,ifelse(po_oldsurv$variable == "survive9", 2015,ifelse(po_oldsurv$variable  == "survive10", 2016, NA)))))))))
# View(po_oldsurv)

po_oldgrow <- POSY_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("TotTillers1", "TotTillers3","TotTillers4", 
                       "TotTillers5", "TotTillers6","TotTillers7", 
                       "TotTillers8", "TotTillers9", "TotTillers10"), 
       value.name = "size") 
po_oldgrow$year<- ifelse(po_oldgrow$variable == "TotTillers1", 2008, ifelse(po_oldgrow$variable  == "TotTillers3", 2009, ifelse(po_oldgrow$variable  == "TotTillers4", 2010, ifelse(po_oldgrow$variable  == "TotTillers5", 2011, ifelse(po_oldgrow$variable  == "TotTillers6", 2012, ifelse(po_oldgrow$variable  == "TotTillers7", 2013, ifelse(po_oldgrow$variable == "TotTillers8", 2014, ifelse(po_oldgrow$variable == "TotTillers9", 2015, ifelse(po_oldgrow$variable  == "TotTillers10", 2016, NA)))))))))
# View(po_oldgrow)

po_oldflw <- POSY_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n",
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("FlwTillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
po_oldflw$year<- ifelse(po_oldflw$variable == "FlwTillers1", 2008, ifelse(po_oldflw$variable  == "FlwTillers3", 2009, ifelse(po_oldflw$variable  == "FlwTillers4", 2010, ifelse(po_oldflw$variable  == "FlwTillers5", 2011, ifelse(po_oldflw$variable  == "FlwTillers6", 2012, ifelse(po_oldflw$variable  == "FlwTillers7", 2013,ifelse(po_oldflw$variable == "FlwTillers8", 2014,ifelse(po_oldflw$variable == "FlwTillers9", 2015,ifelse(po_oldflw$variable  == "FlwTillers10", 2016, NA)))))))))
# View(po_oldflw)


po_oldmerge_sg <- merge(po_oldsurv, po_oldgrow, by = c( "plot","pos", "tag", "Endo", 
                                                  "Loc'n", "Birth Year", "TRT",
                                                  "Plant", "year"))
# View(po_oldmerge_sg)

po_oldmerge_sgf <- merge(po_oldmerge_sg, po_oldflw, by = c( "plot","pos", "tag", "Endo", 
                                                      "Loc'n", "Birth Year", "TRT",
                                                      "Plant", "year"))
# View(po_ldmerge_sgf)

# getting a dataframe with t and t_1
po_oldmerge_t1 <-po_oldmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(po_oldmerge_t1)

po_oldmerge_t <-po_oldmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(po_oldmerge_t)

po_oldmerge <- po_oldmerge_t1 %>% 
  full_join(po_oldmerge_t, by = c("plot","pos", "tag", "Endo", 
                                "Loc'n", "Birth Year", "TRT",
                                "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(po_oldmerge)


# Combining measurements across years for the “Old” POSY recruits data --------

## Combining data for recruits across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
po_roldsurv <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive09 = "NA") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("survive09", "survive10", "survive11", "survive12", "survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
po_roldsurv$year<- ifelse(po_roldsurv$variable == "survive09", 2009, ifelse(po_roldsurv$variable == "survive10", 2010, ifelse(po_roldsurv$variable == "survive11", 2011, ifelse(po_roldsurv$variable  == "survive12", 2012, ifelse(po_roldsurv$variable  == "survive13", 2013, ifelse(po_roldsurv$variable  == "Survive14", 2014, ifelse(po_roldsurv$variable  == "Survive15", 2015, ifelse(po_roldsurv$variable  == "Survive16", 2016, NA))))))))
# View(po_roldsurv)


po_roldgrow <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller09", "TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
po_roldgrow$year<- ifelse(po_roldgrow$variable == "TOTtiller09", 2009, ifelse(po_roldgrow$variable == "TOTtiller10", 2010, ifelse(po_roldgrow$variable == "TOTtiller11", 2011, ifelse(po_roldgrow$variable  == "TOTtiller12", 2012, ifelse(po_roldgrow$variable  == "TOTtiller13", 2013, ifelse(po_roldgrow$variable  == "TOTtiller14", 2014, ifelse(po_roldgrow$variable  == "TOTtiller15", 2015, ifelse(po_roldgrow$variable  == "TOTtiller16", 2016, NA))))))))
# View(po_roldgrow)

po_roldflw <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(FLWtiller09 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller09", "FLWtiller10","FLWtiller11", "FLWTiller12", 
                       "FLWTiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_roldflw$year<- ifelse(po_roldflw$variable == "FLWtiller09", 2009, ifelse(po_roldflw$variable == "FLWtiller10", 2010, ifelse(po_roldflw$variable == "FLWtiller11", 2011, ifelse(po_roldflw$variable  == "FLWTiller12", 2012, ifelse(po_roldflw$variable  == "FLWTiller13", 2013, ifelse(po_roldflw$variable  == "FLWtiller14", 2014, ifelse(po_roldflw$variable  == "FLWtiller15", 2015, ifelse(po_roldflw$variable  == "FLWtiller16", 2016, NA))))))))
# View(po_roldflw)

po_roldmerge_sg <- merge(po_roldsurv, po_roldgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_roldmerge_sg)

po_roldmerge_sgf <- merge(po_roldmerge_sg, po_roldflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_roldmerge_sgf)

## getting a dataframe with time t and t_1
po_roldmerge_t1 <-po_roldmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(po_roldmerge_t1)

po_roldmerge_t <-po_roldmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(po_roldmerge_t)

po_roldmerge <- po_roldmerge_t1 %>% 
  full_join(po_roldmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(po_roldmerge)






# Combining the old and new and original and recruit POSY dataframes ---------
po_merge <- po_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                   "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                   "year_t", "size_t", "flw_t")]
po_ldmerge <- po_oldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "size_t", "flw_t")]
po_rmerge <- po_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                   "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                   "year_t", "size_t", "flw_t")]
po_roldmerge <- po_roldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "size_t", "flw_t")]

POSY <- po_merge %>% 
  rbind(po_oldmerge) %>% 
  rbind(po_roldmerge) %>% 
  rbind(po_rmerge) %>% 
  mutate(species = "POSY") %>% 
  mutate(surv_t1 = ifelse(surv_t1 == 'NA', NA, surv_t1))
POSY <- POSY[!(is.na(POSY$surv_t1)),]

View(POSY)






# Combining measurements across years for the LOAR data ------------------

## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

lsurv <- LOAR_data %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("survive1", "survive2", "survive3","survive5",
                       "survive6(2012)", "survive7(2013)", "survive8(2014)", 
                       "survive9(2015)", "survive10(2016)"),
       value.name = "surv") 
lsurv$year<- ifelse(lsurv$variable == "survive1", 2008, ifelse(lsurv$variable  == "survive2", 2009, ifelse(lsurv$variable  == "survive3", 2010, ifelse(lsurv$variable  == "survive5", 2011, ifelse(lsurv$variable  == "survive6(2012)", 2012, ifelse(lsurv$variable  == "survive7(2013)", 2013,ifelse(lsurv$variable == "survive8(2014)", 2014,ifelse(lsurv$variable == "survive9(2015)", 2015,ifelse(lsurv$variable  == "survive10(2016)", 2016, NA)))))))))
# View(lsurv)

lgrow <- LOAR_data %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("TOTtillers1", "TOTtiller2", "TOTtillers3",
                      "TOTtillers5", "TOTtillers6","TOTtillers7", 
                      "TOTtillers8", "TOTtillers9", "TOTtillers10"), 
       value.name = "size") 
lgrow$year<- ifelse(lgrow$variable == "TOTtillers1", 2008, ifelse(lgrow$variable  == "TOTtiller2", 2009, ifelse(lgrow$variable  == "TOTtillers3", 2010, ifelse(lgrow$variable  == "TOTtillers5", 2011, ifelse(lgrow$variable  == "TOTtillers6", 2012, ifelse(lgrow$variable  == "TOTtillers7", 2013, ifelse(lgrow$variable == "TOTtillers8", 2014, ifelse(lgrow$variable == "TOTtillers9", 2015, ifelse(lgrow$variable  == "TOTtillers10", 2016, NA)))))))))
# View(lgrow)

lflw <- LOAR_data %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("FlwTillers1", "FLWTiller2", "FLWTillers3", 
                       "FLWTillers5", "FLWTillers6", "FLWTillers7", 
                       "FLWTillers8", "FLWTillers9", "FLWTillers10"), 
       value.name = "flw") 
lflw$year<- ifelse(lflw$variable == "FlwTillers1", 2008, ifelse(lflw$variable  == "FLWTiller2", 2009, ifelse(lflw$variable  == "FLWTillers3", 2010, ifelse(lflw$variable  == "FLWTillers5", 2011, ifelse(lflw$variable  == "FLWTillers6", 2012, ifelse(lflw$variable  == "FLWTillers7", 2013,ifelse(lflw$variable == "FLWTillers8", 2014,ifelse(lflw$variable == "FLWTillers9", 2015,ifelse(lflw$variable  == "FLWTillers10", 2016, NA)))))))))
# View(lflw)

l_merge_sg <- merge(lsurv, lgrow, by = c( "plot", "pos", "tag", "Endo", 
                                         "Loc'n", "Birth Year", "TRT",
                                         "Plant", "year"))
# View(l_merge_sg)

l_merge_sgf <- merge(l_merge_sg, lflw, by = c( "plot", "pos", "tag", "Endo", 
                                             "Loc'n", "Birth Year", "TRT",
                                             "Plant", "year"))
# View(l_merge_sgf)

# getting a dataframe with t and t_1
l_merge_t1 <-l_merge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(l_merge_t1)

l_merge_t <-l_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(l_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
l_merge <- l_merge_t1 %>% 
  full_join(l_merge_t, by = c("plot", "pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(l_merge)


# Combining measurements across years for the LOAR recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
l_rsurv <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  mutate(survive10 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
l_rsurv$year<- ifelse(l_rsurv$variable == "survive10", 2010, ifelse(l_rsurv$variable == "Survive11", 2011, ifelse(l_rsurv$variable  == "Survive12", 2012, ifelse(l_rsurv$variable  == "Survive13", 2013, ifelse(l_rsurv$variable  == "Survive14", 2014, ifelse(l_rsurv$variable  == "Survive15", 2015, ifelse(l_rsurv$variable  == "Survive16", 2016, NA)))))))
# View(l_rsurv)

l_rgrow <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TotTillers10","Tot11", "Tot12", 
                       "Tot13", "Tot14", "Tot15", 
                       "Tot16"),
       value.name = "size") 
l_rgrow$year<- ifelse(l_rgrow$variable == "TotTillers10", 2010, ifelse(l_rgrow$variable == "Tot11", 2011, ifelse(l_rgrow$variable  == "Tot12", 2012, ifelse(l_rgrow$variable  == "Tot13", 2013, ifelse(l_rgrow$variable  == "Tot14", 2014, ifelse(l_rgrow$variable  == "Tot15", 2015, ifelse(l_rgrow$variable  == "Tot16", 2016, NA)))))))
# View(l_rgrow)

l_rflw <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FlwTillers10","Flw11", "Flw12", 
                       "Flw13", "Flw14", "Flw15",
                       "Flw16"),
       value.name = "flw") 
l_rflw$year<- ifelse(l_rflw$variable == "FlwTillers10", 2010, ifelse(l_rflw$variable == "Flw11", 2011, ifelse(l_rflw$variable  == "Flw12", 2012, ifelse(l_rflw$variable  == "Flw13", 2013, ifelse(l_rflw$variable  == "Flw14", 2014, ifelse(l_rflw$variable  == "Flw15", 2015, ifelse(l_rflw$variable  == "Flw16", 2016, NA)))))))
# View(l_rflw)

l_rmerge_sg <- merge(l_rsurv, l_rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(l_rmerge_sg)

l_rmerge_sgf <- merge(l_rmerge_sg, l_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(l_rmerge_sgf)

## getting a dataframe with time t and t_1
l_rmerge_t1 <-l_rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(l_rmerge_t1)

l_rmerge_t <-l_rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(l_rmerge_t)

l_rmerge <- l_rmerge_t1 %>% 
  full_join(l_rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(l_rmerge)






# Combining the  original and recruit LOAR dataframes ---------
l_merge <- l_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                   "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                   "year_t", "size_t", "flw_t")]

l_rmerge <- l_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                   "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                   "year_t", "size_t", "flw_t")]


LOAR <- l_merge %>% 
  rbind(l_rmerge) %>% 
  mutate(species = "LOAR")
LOAR <- LOAR[!(is.na(LOAR$surv_t1)),]

View(LOAR)












# Combining measurements across years for the FESU data ------------------

## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

fsurv <- FESU_data %>%
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("survive1_2008", "Survive2", "Survive3","survive4",
                       "survive5", "survive6", "survive7", 
                       "survive8", "survive9"),
       value.name = "surv") 
fsurv$year<- ifelse(fsurv$variable == "survive1_2008", 2008, ifelse(fsurv$variable  == "Survive2", 2009, ifelse(fsurv$variable  == "Survive3", 2010, ifelse(fsurv$variable  == "survive4", 2011, ifelse(fsurv$variable  == "survive5", 2012, ifelse(fsurv$variable  == "survive6", 2013,ifelse(fsurv$variable == "survive7", 2014,ifelse(fsurv$variable == "survive8", 2015,ifelse(fsurv$variable  == "survive9", 2016, NA)))))))))
# View(fsurv)

fgrow <- FESU_data %>% 
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("TotTillers1", "TotTillers2", "TotTillers3",
                       "TotTillers4", "TotTillers5","TotTillers6", 
                       "TotTillers7", "TotTillers8", "TotTillers9"), 
       value.name = "size") 
fgrow$year<- ifelse(fgrow$variable == "TotTillers1", 2008, ifelse(fgrow$variable  == "TotTillers2", 2009, ifelse(fgrow$variable  == "TotTillers3", 2010, ifelse(fgrow$variable  == "TotTillers4", 2011, ifelse(fgrow$variable  == "TotTillers5", 2012, ifelse(fgrow$variable  == "TotTillers6", 2013, ifelse(fgrow$variable == "TotTillers7", 2014, ifelse(fgrow$variable == "TotTillers8", 2015, ifelse(fgrow$variable  == "TotTillers9", 2016, NA)))))))))
# View(fgrow)

fflw <- FESU_data %>% 
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(FlwTillers1 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("FlwTillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
fflw$year<- ifelse(fflw$variable == "FlwTillers1", 2008, ifelse(fflw$variable  == "FlwTillers2", 2009, ifelse(fflw$variable  == "FlwTillers3", 2010, ifelse(fflw$variable  == "FlwTillers4", 2011, ifelse(fflw$variable  == "FlwTillers5", 2012, ifelse(fflw$variable  == "FlwTillers6", 2013,ifelse(fflw$variable == "FlwTillers7", 2014,ifelse(fflw$variable == "FlwTillers8", 2015,ifelse(fflw$variable  == "FlwTillers9", 2016, NA)))))))))
# View(fflw)

f_merge_sg <- merge(fsurv, fgrow, by = c( "plot", "pos", "tag", "Endo", 
                                          "Loc'n", "Birth Year", "TRT",
                                          "Plant", "year"))
# View(f_merge_sg)

f_merge_sgf <- merge(f_merge_sg, fflw, by = c( "plot", "pos", "tag", "Endo", 
                                               "Loc'n", "Birth Year", "TRT",
                                               "Plant", "year"))
# View(f_merge_sgf)

# getting a dataframe with t and t_1
f_merge_t1 <-f_merge_sgf %>%
  filter(year != min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(f_merge_t1)

f_merge_t <-f_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(f_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
f_merge <- f_merge_t1 %>% 
  full_join(f_merge_t, by = c("plot", "pos", "tag", "Endo", 
                              "Loc'n", "Birth Year", "TRT",
                              "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(f_merge)


# Combining measurements across years for the FESU recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
f_rsurv <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>%
  mutate(survive10 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("survive10", "surv11", "surv12", "Surv13", "Surv14", 
                       "Surv15", "Surv16"),
       value.name = "surv") 
f_rsurv$year<- ifelse(f_rsurv$variable == "survive10", 2010, ifelse(f_rsurv$variable == "surv11", 2011, ifelse(f_rsurv$variable  == "surv12", 2012, ifelse(f_rsurv$variable  == "Surv13", 2013, ifelse(f_rsurv$variable  == "Surv14", 2014, ifelse(f_rsurv$variable  == "Surv15", 2015, ifelse(f_rsurv$variable  == "Surv16", 2016, NA)))))))
# View(f_rsurv)

f_rgrow <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOT10","TOT11", "TOT12", 
                       "TOT13", "TOT14", "TOT15", 
                       "TOT16"),
       value.name = "size") 
f_rgrow$year<- ifelse(f_rgrow$variable == "TOT10", 2010, ifelse(f_rgrow$variable == "TOT11", 2011, ifelse(f_rgrow$variable  == "TOT12", 2012, ifelse(f_rgrow$variable  == "TOT13", 2013, ifelse(f_rgrow$variable  == "TOT14", 2014, ifelse(f_rgrow$variable  == "TOT15", 2015, ifelse(f_rgrow$variable  == "TOT16", 2016, NA)))))))
# View(f_rgrow)

f_rflw <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLW10","FLW11", "FLW12", 
                       "FLW13", "FLW14", "FLW15",
                       "FLW16"),
       value.name = "flw") 
f_rflw$year<- ifelse(f_rflw$variable == "FLW10", 2010, ifelse(f_rflw$variable == "FLW11", 2011, ifelse(f_rflw$variable  == "FLW12", 2012, ifelse(f_rflw$variable  == "FLW13", 2013, ifelse(f_rflw$variable  == "FLW14", 2014, ifelse(f_rflw$variable  == "FLW15", 2015, ifelse(f_rflw$variable  == "FLW16", 2016, NA)))))))
# View(f_rflw)

f_rmerge_sg <- merge(f_rsurv, f_rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(f_rmerge_sg)

f_rmerge_sgf <- merge(f_rmerge_sg, f_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(f_rmerge_sgf)

## getting a dataframe with time t and t_1
f_rmerge_t1 <-f_rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(f_rmerge_t1)

f_rmerge_t <-f_rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(f_rmerge_t)

f_rmerge <- f_rmerge_t1 %>% 
  full_join(f_rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
 # View(f_rmerge)











# Combining the original and recruit FESU dataframes ----------------------
f_merge <- f_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                     "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                     "year_t", "size_t", "flw_t")]

f_rmerge <- f_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                       "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                       "year_t", "size_t", "flw_t")]

FESU <- f_merge %>% 
  rbind(f_rmerge) %>% 
  mutate(species = "FESU")
FESU <- FESU[!(is.na(FESU$surv_t1)),]

View(FESU)

# Combining measurements across years for the ELVI data ------------------

## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

elvisurv <- ELVI_data %>%
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("Survive08", "Survive09", "Survive10","Survive11",
                       "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
elvisurv$year<- ifelse(elvisurv$variable == "Survive08", 2008, ifelse(elvisurv$variable  == "Survive09", 2009, ifelse(elvisurv$variable  == "Survive10", 2010, ifelse(elvisurv$variable  == "Survive11", 2011, ifelse(elvisurv$variable  == "Survive12", 2012, ifelse(elvisurv$variable  == "Survive13", 2013,ifelse(elvisurv$variable == "Survive14", 2014,ifelse(elvisurv$variable == "Survive15", 2015,ifelse(elvisurv$variable  == "Survive16", 2016, NA)))))))))
# View(elvisurv)

elvigrow <- ELVI_data %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("TotTillers07", "TotTillers08", "TotTillers09", "Tottillers10",
                       "TotTillers11", "TotTillers12","TotTillers13", 
                       "TotTillers14", "TotTillers15", "TotTillers16"), 
       value.name = "size") 
elvigrow$year<- ifelse(elvigrow$variable == "TotTillers07", 2007, ifelse(elvigrow$variable == "TotTillers08", 2008, ifelse(elvigrow$variable  == "TotTillers09", 2009, ifelse(elvigrow$variable  == "Tottillers10", 2010, ifelse(elvigrow$variable  == "TotTillers11", 2011, ifelse(elvigrow$variable  == "TotTillers12", 2012, ifelse(elvigrow$variable  == "TotTillers13", 2013, ifelse(elvigrow$variable == "TotTillers14", 2014, ifelse(elvigrow$variable == "TotTillers15", 2015, ifelse(elvigrow$variable  == "TotTillers16", 2016, NA))))))))))
# View(elvigrow)

elviflw <- ELVI_data %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(FlwTillers1 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("FlwTillers08", "FlwTillers09", "FlwTillers10", 
                       "FlwTillers11", "FlwTillers12", "FlwTillers13", 
                       "FlwTillers14", "FlwTillers15", "FlwTillers16"), 
       value.name = "flw") 
elviflw$year<- ifelse(elviflw$variable == "FlwTillers08", 2008, ifelse(elviflw$variable  == "FlwTillers09", 2009, ifelse(elviflw$variable  == "FlwTillers10", 2010, ifelse(elviflw$variable  == "FlwTillers11", 2011, ifelse(elviflw$variable  == "FlwTillers12", 2012, ifelse(elviflw$variable  == "FlwTillers13", 2013,ifelse(elviflw$variable == "FlwTillers14", 2014,ifelse(elviflw$variable == "FlwTillers15", 2015,ifelse(elviflw$variable  == "FlwTillers16", 2016, NA)))))))))
# View(elviflw)

elvi_merge_sg <- merge(elvisurv, elvigrow, by = c( "plot", "pos", "tag", "Endo", 
                                          "Loc'n", "Birth Year", "TRT",
                                          "Plant", "year"))
# View(elvi_merge_sg)

elvi_merge_sgf <- merge(elvi_merge_sg, elviflw, by = c( "plot", "pos", "tag", "Endo", 
                                               "Loc'n", "Birth Year", "TRT",
                                               "Plant", "year"))
# View(elvi_merge_sgf)

# getting a dataframe with t and t_1
elvi_merge_t1 <-elvi_merge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(elvi_merge_t1)

elvi_merge_t <-elvi_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(elvi_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
elvi_merge <- elvi_merge_t1 %>% 
  full_join(elvi_merge_t, by = c("plot", "pos", "tag", "Endo", 
                              "Loc'n", "Birth Year", "TRT",
                              "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
 # View(elvi_merge)


# Combining measurements across years for the ELVI recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
elvi_rsurv <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("Survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15"),
       value.name = "surv") 
elvi_rsurv$year<- ifelse(elvi_rsurv$variable == "Survive10", 2010, ifelse(elvi_rsurv$variable == "Survive11", 2011, ifelse(elvi_rsurv$variable  == "Survive12", 2012, ifelse(elvi_rsurv$variable  == "Survive13", 2013, ifelse(elvi_rsurv$variable  == "Survive14", 2014, ifelse(elvi_rsurv$variable  == "Survive15", 2015, NA))))))
# View(elvi_rsurv)

elvi_rgrow <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TotTillers10","TotTillers11", "TotTillers12", 
                       "TotTillers13", "TotTillers14", "TotTillers15"),
       value.name = "size") 
elvi_rgrow$year<- ifelse(elvi_rgrow$variable == "TotTillers10", 2010, ifelse(elvi_rgrow$variable == "TotTillers11", 2011, ifelse(elvi_rgrow$variable  == "TotTillers12", 2012, ifelse(elvi_rgrow$variable  == "TotTillers13", 2013, ifelse(elvi_rgrow$variable  == "TotTillers14", 2014, ifelse(elvi_rgrow$variable  == "TotTillers15", 2015, NA))))))
# View(elvi_rgrow)

elvi_rflw <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15"),
       value.name = "flw") 
elvi_rflw$year<- ifelse(elvi_rflw$variable == "FlwTillers10", 2010, ifelse(elvi_rflw$variable == "FlwTillers11", 2011, ifelse(elvi_rflw$variable  == "FlwTillers12", 2012, ifelse(elvi_rflw$variable  == "FlwTillers13", 2013, ifelse(elvi_rflw$variable  == "FlwTillers14", 2014, ifelse(elvi_rflw$variable  == "FlwTillers15", 2015, ifelse(elvi_rflw$variable  == "FlwTillers16", 2016, NA)))))))
# View(elvi_rflw)

elvi_rmerge_sg <- merge(elvi_rsurv, elvi_rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(elvi_rmerge_sg)

elvi_rmerge_sgf <- merge(elvi_rmerge_sg, elvi_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(elvi_rmerge_sgf)

## getting a dataframe with time t and t_1
elvi_rmerge_t1 <-elvi_rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(elvi_rmerge_t1)

elvi_rmerge_t <-elvi_rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(elvi_rmerge_t)

elvi_rmerge <- elvi_rmerge_t1 %>% 
  full_join(elvi_rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(elvi_rmerge)










# Combining the  original and recruit ELVI dataframes ---------
elvi_merge <- elvi_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                     "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                     "year_t", "size_t", "flw_t")]

elvi_rmerge <- elvi_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                       "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                       "year_t", "size_t", "flw_t")]


ELVI <- elvi_merge %>% 
  rbind(elvi_rmerge) %>% 
  mutate(species = "ELVI")
ELVI <- ELVI[!(is.na(ELVI$surv_t1)),]

View(ELVI)















# Combining measurements across years for the ELRI data ------------------

## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

elrisurv <- ELRI_data %>%
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>% 
  mutate(Survive07 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("Survive07", "survive08", "Survive09", "Survive10","Survive11",
                       "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
elrisurv$year<- ifelse(elrisurv$variable == "Survive07", 2007, ifelse(elrisurv$variable == "survive08", 2008, ifelse(elrisurv$variable  == "Survive09", 2009, ifelse(elrisurv$variable  == "Survive10", 2010, ifelse(elrisurv$variable  == "Survive11", 2011, ifelse(elrisurv$variable  == "Survive12", 2012, ifelse(elrisurv$variable  == "Survive13", 2013,ifelse(elrisurv$variable == "Survive14", 2014,ifelse(elrisurv$variable == "Survive15", 2015,ifelse(elrisurv$variable  == "Survive16", 2016, NA))))))))))
# View(elrisurv)

elrigrow <- ELRI_data %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("TotTillers07", "TotTillers08", "TotTillers09", "TotTillers10",
                       "TotTiller11", "TotTiller12","TotTiller13", 
                       "TotTiller14", "TotTiller15", "TotTillers16"), 
       value.name = "size") 
elrigrow$year<- ifelse(elrigrow$variable == "TotTillers07", 2007, ifelse(elrigrow$variable == "TotTillers08", 2008, ifelse(elrigrow$variable  == "TotTillers09", 2009, ifelse(elrigrow$variable  == "TotTillers10", 2010, ifelse(elrigrow$variable  == "TotTiller11", 2011, ifelse(elrigrow$variable  == "TotTiller12", 2012, ifelse(elrigrow$variable  == "TotTiller13", 2013, ifelse(elrigrow$variable == "TotTiller14", 2014, ifelse(elrigrow$variable == "TotTiller15", 2015, ifelse(elrigrow$variable  == "TotTillers16", 2016, NA))))))))))
# View(elrigrow)

elriflw <- ELRI_data %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  mutate(FlwTillers1 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("FlwTillers08", "FlwTillers09", "FLwTillers10", 
                       "FlwTiller11", "FlwTiller12", "FlwTiller13", 
                       "FlwTiller14", "FlwTiller15", "FlwTillers16"), 
       value.name = "flw") 
elriflw$year<- ifelse(elriflw$variable == "FlwTillers08", 2008, ifelse(elriflw$variable  == "FlwTillers09", 2009, ifelse(elriflw$variable  == "FLwTillers10", 2010, ifelse(elriflw$variable  == "FlwTiller11", 2011, ifelse(elriflw$variable  == "FlwTiller12", 2012, ifelse(elriflw$variable  == "FlwTiller13", 2013,ifelse(elriflw$variable == "FlwTiller14", 2014,ifelse(elriflw$variable == "FlwTiller15", 2015,ifelse(elriflw$variable  == "FlwTillers16", 2016, NA)))))))))
# View(elriflw)

elri_merge_sg <- merge(elrisurv, elrigrow, by = c( "plot", "pos", "tag", "Endo", 
                                                   "Loc'n", "Birth Year", "TRT",
                                                   "Plant", "year"))
# View(elri_merge_sg)

elri_merge_sgf <- merge(elri_merge_sg, elriflw, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Loc'n", "Birth Year", "TRT",
                                                        "Plant", "year"))
# View(elri_merge_sgf)

# getting a dataframe with t and t_1
elri_merge_t1 <-elri_merge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(elri_merge_t1)

elri_merge_t <-elri_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(elri_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
elri_merge <- elri_merge_t1 %>% 
  full_join(elri_merge_t, by = c("plot", "pos", "tag", "Endo", 
                                 "Loc'n", "Birth Year", "TRT",
                                 "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(elri_merge)


# Combining measurements across years for the ELRI recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
elri_rsurv <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("surv10", "surv11", "surv12", "Surv13", "Surv14", 
                       "Surv15", "Survive16"),
       value.name = "surv") 
elri_rsurv$year<- ifelse(elri_rsurv$variable == "surv10", 2010, ifelse(elri_rsurv$variable == "surv11", 2011, ifelse(elri_rsurv$variable  == "surv12", 2012, ifelse(elri_rsurv$variable  == "Surv13", 2013, ifelse(elri_rsurv$variable  == "Surv14", 2014, ifelse(elri_rsurv$variable  == "Surv15", 2015, ifelse(elri_rsurv$variable == "Survive16", 2016, NA)))))))
# View(elri_rsurv)

elri_rgrow <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller09", "TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOT13", "TOT14", "TOT15", "TotTillers16"),
       value.name = "size") 
elri_rgrow$year<- ifelse(elri_rgrow$variable == "TOTtiller09", 2009, ifelse(elri_rgrow$variable == "TOTtiller10", 2010, ifelse(elri_rgrow$variable == "TOTtiller11", 2011, ifelse(elri_rgrow$variable  == "TOTtiller12", 2012, ifelse(elri_rgrow$variable  == "TOT13", 2013, ifelse(elri_rgrow$variable  == "TOT14", 2014, ifelse(elri_rgrow$variable  == "TOT15", 2015, ifelse(elri_rgrow$variable == "TotTillers16", 2016, NA))))))))
# View(elri_rgrow)

elri_rflw <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLW13", "FLW14", "FLW15", "FlwTillers16"),
       value.name = "flw") 
elri_rflw$year<- ifelse(elri_rflw$variable == "FLWtiller09", 2009, ifelse(elri_rflw$variable == "FLWtiller10", 2010, ifelse(elri_rflw$variable == "FLWtiller11", 2011, ifelse(elri_rflw$variable  == "FLWtiller12", 2012, ifelse(elri_rflw$variable  == "FLW13", 2013, ifelse(elri_rflw$variable  == "FLW14", 2014, ifelse(elri_rflw$variable  == "FLW15", 2015, ifelse(elri_rflw$variable  == "FlwTillers16", 2016, NA))))))))
# View(elri_rflw)

elri_rmerge_sg <- merge(elri_rsurv, elri_rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(elri_rmerge_sg)

elri_rmerge_sgf <- merge(elri_rmerge_sg, elri_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(elri_rmerge_sgf)

## getting a dataframe with time t and t_1
elri_rmerge_t1 <-elri_rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(elri_rmerge_t1)

elri_rmerge_t <-elri_rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(elri_rmerge_t)

elri_rmerge <- elri_rmerge_t1 %>% 
  full_join(elri_rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(elri_rmerge)










# Combining the  original and recruit ELRI dataframes ---------
elri_merge <- elri_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                           "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                           "year_t", "size_t", "flw_t")]

elri_rmerge <- elri_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                             "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                             "year_t", "size_t", "flw_t")]


ELRI <- elri_merge %>% 
  rbind(elri_rmerge) %>% 
  mutate(species = "ELRI")
ELRI <- ELRI[!(is.na(ELRI$surv_t1)),]

View(ELRI)



















# Combining measurements across years for the AGPE data ------------------

## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

agpesurv <- AGPE_data %>%
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate(Survive07 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("survive1", "survive2", "survive3","survive4",
                       "survive5", "Survive6", "Survive7", 
                       "Survive8", "Survive9"),
       value.name = "surv") 
agpesurv$year<- ifelse(agpesurv$variable == "survive1", 2008, ifelse(agpesurv$variable  == "survive2", 2009, ifelse(agpesurv$variable  == "survive3", 2010, ifelse(agpesurv$variable  == "survive4", 2011, ifelse(agpesurv$variable  == "survive5", 2012, ifelse(agpesurv$variable  == "Survive6", 2013,ifelse(agpesurv$variable == "Survive7", 2014,ifelse(agpesurv$variable == "Survive8", 2015,ifelse(agpesurv$variable  == "Survive9", 2016, NA)))))))))
# View(agpesurv)

agpegrow <- AGPE_data %>% 
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("ToTillers07", "Tottillers1", "TotTillers2", 
                       "TOTtiller3", "TOTtiller4", "TotTiller5", 
                       "TotTillers6", "TotTillers7", "TotTillers8", 
                       "TotTillers9"), 
       value.name = "size") 
agpegrow$year<- ifelse(agpegrow$variable == "ToTillers07", 2007, ifelse(agpegrow$variable == "TotTillers1", 2008, ifelse(agpegrow$variable  == "TotTillers2", 2009, ifelse(agpegrow$variable  == "TOTtiller3", 2010, ifelse(agpegrow$variable  == "TOTtiller4", 2011, ifelse(agpegrow$variable  == "TotTiller5", 2012, ifelse(agpegrow$variable  == "TotTillers6", 2013, ifelse(agpegrow$variable == "TotTillers7", 2014, ifelse(agpegrow$variable == "TotTillers8", 2015, ifelse(agpegrow$variable  == "TotTillers9", 2016, NA))))))))))
# View(agpegrow)

agpeflw <- AGPE_data %>% 
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(FlwTillers1 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("FlwTillers1", "FlwTillers2", "FLWTiller3", 
                       "FLWTiller4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
agpeflw$year<- ifelse(agpeflw$variable == "FlwTillers1", 2008, ifelse(agpeflw$variable  == "FlwTillers2", 2009, ifelse(agpeflw$variable  == "FLWTiller3", 2010, ifelse(agpeflw$variable  == "FLWtiller4", 2011, ifelse(agpeflw$variable  == "FlwTillers5", 2012, ifelse(agpeflw$variable  == "FlwTillers6", 2013,ifelse(agpeflw$variable == "FlwTillers7", 2014,ifelse(agpeflw$variable == "FlwTillers8", 2015,ifelse(agpeflw$variable  == "FlwTillers9", 2016, NA)))))))))
# View(agpeflw)

agpe_merge_sg <- merge(agpesurv, agpegrow, by = c( "plot", "pos", "tag", "Endo", 
                                                   "Loc'n", "Birth Year", "TRT",
                                                   "Plant", "year"))
# View(agpe_merge_sg)

agpe_merge_sgf <- merge(agpe_merge_sg, agpeflw, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Loc'n", "Birth Year", "TRT",
                                                        "Plant", "year"))
View(agpe_merge_sgf)

# getting a dataframe with t and t_1
agpe_merge_t1 <-agpe_merge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(agpe_merge_t1)

agpe_merge_t <-agpe_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(agpe_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
agpe_merge <- agpe_merge_t1 %>% 
  full_join(agpe_merge_t, by = c("plot", "pos", "tag", "Endo", 
                                 "Loc'n", "Birth Year", "TRT",
                                 "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(agpe_merge)


# Combining measurements across years for the AGPE recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
agpe_rsurv <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("Survive09", "Survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
agpe_rsurv$year<- ifelse(agpe_rsurv$variable == "Survive09", 2009, ifelse(agpe_rsurv$variable == "Survive10", 2010, ifelse(agpe_rsurv$variable == "Survive11", 2011, ifelse(agpe_rsurv$variable  == "Survive12", 2012, ifelse(agpe_rsurv$variable  == "Survive13", 2013, ifelse(agpe_rsurv$variable  == "Survive14", 2014, ifelse(agpe_rsurv$variable  == "Survive15", 2015, ifelse(agpe_rsurv$variable == "Survive16", 2016, NA))))))))
# View(agpe_rsurv)

agpe_rgrow <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TotTillers09", "TotTillers10","TotTillers11", "TotTiller12", 
                       "TotTillers13", "TotTillers14", "TotTillers15", "TotTillers16"),
       value.name = "size") 
agpe_rgrow$year<- ifelse(agpe_rgrow$variable == "TotTillers09", 2009, ifelse(agpe_rgrow$variable == "TotTillers10", 2010, ifelse(agpe_rgrow$variable == "TotTillers11", 2011, ifelse(agpe_rgrow$variable  == "TotTiller12", 2012, ifelse(agpe_rgrow$variable  == "TotTillers13", 2013, ifelse(agpe_rgrow$variable  == "TotTillers14", 2014, ifelse(agpe_rgrow$variable  == "TotTillers15", 2015, ifelse(agpe_rgrow$variable == "TotTillers16", 2016, NA))))))))
# View(agpe_rgrow)

agpe_rflw <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FlwTillers09", "FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15", "FlwTillers16"),
       value.name = "flw") 
agpe_rflw$year<- ifelse(agpe_rflw$variable == "FlwTillers09", 2009, ifelse(agpe_rflw$variable == "FlwTillers10", 2010, ifelse(agpe_rflw$variable == "FlwTillers11", 2011, ifelse(agpe_rflw$variable  == "FlwTillers12", 2012, ifelse(agpe_rflw$variable  == "FlwTillers13", 2013, ifelse(agpe_rflw$variable  == "FlwTillers14", 2014, ifelse(agpe_rflw$variable  == "FlwTillers15", 2015, ifelse(agpe_rflw$variable  == "FlwTillers16", 2016, NA))))))))
# View(agpe_rflw)

agpe_rmerge_sg <- merge(agpe_rsurv, agpe_rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(agpe_rmerge_sg)

agpe_rmerge_sgf <- merge(agpe_rmerge_sg, agpe_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(agpe_rmerge_sgf)

## getting a dataframe with time t and t_1
agpe_rmerge_t1 <-agpe_rmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(agpe_rmerge_t1)

agpe_rmerge_t <-agpe_rmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(agpe_rmerge_t)

agpe_rmerge <- agpe_rmerge_t1 %>% 
  full_join(agpe_rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(agpe_rmerge)










# Combining the  original and recruit AGPE dataframes ---------
agpe_merge <- agpe_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                           "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                           "year_t", "size_t", "flw_t")]

agpe_rmerge <- agpe_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                             "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                             "year_t", "size_t", "flw_t")]


AGPE <- agpe_merge %>% 
  rbind(agpe_rmerge) %>% 
  mutate(species = "AGPE")
AGPE <- AGPE[!(is.na(AGPE$surv_t1)),]
View(AGPE)



# Combining all of the species dataframes, without the 2017 data which is in endo_demog_long ----------
# There are still some oddities in this set. For example: The surv_t1 column has 0,1,Na, 1? 0? and XX as entries
LTREB_endodemog <- AGPE %>% 
  rbind(LOAR) %>% 
  rbind(ELVI) %>% 
  rbind(ELRI) %>% 
  rbind(FESU) %>% 
  rbind(POAL) %>%  
  rbind(POSY)
 # View(LTREB_endodemog)




# Pulling out the 2017 data from endo_demog_long --------------------------

endo_demog_long<- read.csv("/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
# View(endo_demog_long)

endo2017 <- endo_demog_long %>% 
  rename("tag" = "id", "Endo" = "endo", "Loc'n" = "quad", 
         "Birth Year" = "birth" ) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA) %>% 
  select("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
         "year_t", "size_t", "species") %>% 
  filter(year_t1 == "2018" | year_t1 == "2017") %>% 
  filter(!is.na(surv_t1))
unique(endo2017$Endo)
endo2017$Endo <- as.integer(ifelse(endo2017$Endo == "1", "1", ifelse(endo2017$Endo == "0", "0", ifelse(endo2017$Endo == "plus", "1", ifelse(endo2017$Endo == "minus", "0", NA)))))
endo2017$origin <- ifelse(endo2017$origin == "O", 0, ifelse(endo2017$origin == "R", 1, NA))
  
endo2017<- endo2017[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                      "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                      "year_t", "size_t", "species")]
str(endo2017)


# Combining the 2017 data with the multi species data frame.

LTREB_endodemog <- LTREB_endodemog %>% 
  select("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
         "year_t", "size_t", "species") %>% 
  rbind(endo2017)
 # View(LTREB_endodemog)

# figuring out an issue with NAs in the size_t
# I am working on how to do : for(tag,) if size_t = NA and surv_t1 = 0, remove all but the lowest year

sizetNA <- LTREB_endodemog[which(is.na(LTREB_endodemog$size_t)),] 
unique(sizetNA$species)
unique(sizetNA$surv_t1)
unique(sizetNA$year_t)
# Now we can check for funky, out of place values
str(LTREB_endodemog)
unique(LTREB_endodemog$`Birth Year`)
unique(LTREB_endodemog$surv_t1)
unique(LTREB_endodemog$size_t)
unique(LTREB_endodemog$size_t1)
unique(LTREB_endodemog$flw_t1)
unique(LTREB_endodemog$Endo)

View(LTREB_endodemog)

# There are survival values with `?` and `XX`. 
# Reassign the surv_t1 columns to be numeric
xx <- LTREB_endodemog[which(LTREB_endodemog$surv_t1 == "XX"),] 
q1 <- LTREB_endodemog[which(LTREB_endodemog$surv_t1 == "1?"),]
q0 <- LTREB_endodemog[which(LTREB_endodemog$surv_t1 == "0?"),]
LTREB_endodemog[which(LTREB_endodemog$surv_t1 == "XX"),] <- "0"
LTREB_endodemog[which(LTREB_endodemog$surv_t1 == "1?"),] <- "1"
LTREB_endodemog[which(LTREB_endodemog$surv_t1 == "0?"),] <- "0"


LTREB_endodemog$surv_t1 <- as.numeric(as.character(LTREB_endodemog$surv_t1))
str(LTREB_endodemog)

# reassign the Birth Year column to be numeric
LTREB_endodemog[which(LTREB_endodemog$'Birth Year' == "Q2010"),]

LTREB_endodemog$`Birth Year` <- as.numeric(as.character(LTREB_endodemog$`Birth Year`))


# There is a non-integer value for size, 5.5. 
# This is labelled in the data as an avg. between two years because the plant was not found in 2012, but was found in 2011 and 2013. I'm leaving it as is for now.

# Reassign the size columns to be numeric
LTREB_endodemog[which(LTREB_endodemog$size_t1 == 5.5),] 
LTREB_endodemog$size_t1 <- as.numeric(as.character(LTREB_endodemog$size_t1))
LTREB_endodemog$size_t <- as.numeric(as.character(LTREB_endodemog$size_t1))

# There is a value for flower tillers of "check tag".
# Reassign the flw_t1 column to be numeric
LTREB_endodemog[which(LTREB_endodemog$flw_t1 == "check tag"),] <- NA

str(LTREB_endodemog)
View(LTREB_endodemog)

# This can still be further cleaned up, but I am going to take this data frame
# as a csv file for now that I can use for individual species models
write_csv(LTREB_endodemog, "LTREB_endodemog.csv")







# Pulling out the seed production estimates. These are not measured for all plants, and so will go into a separate dataframe------------------------------
# Pulling out the seed production estimates for the "New" POAL data --------
pseed <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2008", "seeds_InflA2", "seeds_InflB2", 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))

pseed$year<- ifelse(pseed$variable == "seed2008", 2008,ifelse(pseed$variable == "seeds_InflA2", 2009, ifelse(pseed$variable  == "seeds_InflB2", 2009, ifelse(pseed$variable  == "seed2010", 2010, ifelse(pseed$variable  == "seed2011", 2011, ifelse(pseed$variable  == "seed2012", 2012, ifelse(pseed$variable  == "seed2013", 2013,ifelse(pseed$variable == "seed2014", 2014,ifelse(pseed$variable == "seed2015", 2015,ifelse(pseed$variable  == "seed2016", 2016, NA))))))))))


 # View(pseed)

pspike <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2008", "spikelets_InflA2", "spikelets_InflB2", "spikelets_inflA3", 
                       "spikelets_inflB3", "spikelets_inflC3", "spikelets_inflA4",
                       "spikelets_inflA5", "spikelets_inflA6", "spikelets_InflB6", 
                       "spikelets_inflC6", "spikelets_inflA7", "spikelets_InflB7", 
                       "spikelets_inflC7","spikelets_inflA8", "spikelets_InflB8", 
                       "spikelets_inflC8","spikelets_inflA9", "spikelets_InflB9", 
                       "spikelets_inflC9"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))
pspike$year<- ifelse(pspike$variable == "spike2008", 2008,ifelse(pspike$variable == "spikelets_InflA2", 2009, ifelse(pspike$variable  == "spikelets_InflB2", 2009, ifelse(pspike$variable  == "spikelets_inflA3", 2010, ifelse(pspike$variable  == "spikelets_inflB3", 2010, ifelse(pspike$variable  == "spikelets_inflB3", 2010, ifelse(pspike$variable  == "spikelets_inflA4", 2011,ifelse(pspike$variable == "spikelets_infA5", 2012,ifelse(pspike$variable == "spikelets_inflA6", 2013, ifelse(pspike$variable == "spikelets_InflB6", 2013, ifelse(pspike$variable == "spikelets_inflC6", 2013, ifelse(pspike$variable  == "spikelets_inflA7", 2014, ifelse(pspike$variable == "spikelets_InflB7", 2014, ifelse(pspike$variable == "spikelets_inflC7", 2014, ifelse(pspike$variable == "spikelets_inflA8", 2015, ifelse(pspike$variable == "spikelets_InflB8", 2015, ifelse(pspike$variable ==   "spikelets_inflC8", 2015, ifelse(pspike$variable == "spikelets_inflA9", 2016, ifelse(pspike$variable == "spikelets_InflB9", 2016, ifelse(pspike$variable == "spikelets_inflC9", 2016, NA))))))))))))))))))))
# View(pspike)

# We already have the FlwTiller data within pflw dataframe

pflw <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
pflw$year<- ifelse(pflw$variable == "Flwtillers1", 2008, ifelse(pflw$variable  == "FlwTillers2", 2009, ifelse(pflw$variable  == "FlwTillers3", 2010, ifelse(pflw$variable  == "FlwTillers4", 2011, ifelse(pflw$variable  == "FlwTillers5", 2012, ifelse(pflw$variable  == "FlwTillers6", 2013,ifelse(pflw$variable == "FlwTillers7", 2014,ifelse(pflw$variable == "FlwTillers8", 2015,ifelse(pflw$variable  == "FlwTillers9", 2016, NA)))))))))

# View(pflw)

pseedmerge_ss <- left_join(pspike, pseed, by = c( "plot", "pos", "tag", "Endo", 
                                         "Birth Year","year", "tillerid"))
# View(pseedmerge_ss)

pseedmerge_ssf <- merge(pseedmerge_ss, pflw, by = c("plot", "pos", "tag", "Endo", 
                                                    "Birth Year","year"), all = TRUE)
# View(pseedmerge_ssf)



# Pulling out the seed, spikelet and flw tiller info from the POAL New Recruits
rseed <-POAL_data_r %>% 
  mutate("Loc'n" = NA, "TRT" = NA) %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2010","seed2011", "seed2012", "seed2013", "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))
rseed$year<- ifelse(rseed$variable == "seed2010", 2010, ifelse(rseed$variable  == "seed2011", 2011, ifelse(rseed$variable  == "seed2012", 2012, ifelse(rseed$variable  == "seed2013", 2013, ifelse(rseed$variable  == "seed2014", 2014, ifelse(rseed$variable  == "seed2015", 2015, ifelse(rseed$variable  == "seed2016", 2016, NA)))))))

# View(rseed)

rspike <- POAL_data_r %>% 
  mutate("Loc'n" = NA, "TRT" = NA) %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2010", "spike2011", "spikelets_inflA5", "spikelets_inflA5__1", "spikelets_infl14", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",))
rspike$year<- ifelse(rspike$variable == "spike2010", 2010, ifelse(rspike$variable == "spike2011", 2011, ifelse(rspike$variable == "spikelets_inflA5", 2012, ifelse(rspike$variable  == "spikelets_inflA5__1", 2013, ifelse(rspike$variable  == "spikelets_infl14", 2014, ifelse(rspike$variable  == "spikelets_infl15", 2015, ifelse(rspike$variable  == "spikelets_infl16", 2016, NA)))))))
# View(rspike)

# We already have the FlwTiller data within rflw dataframe
rflw <- POAL_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
rflw$year<- ifelse(rflw$variable == "FLWtiller10", 2010, ifelse(rflw$variable == "FLWtiller11", 2011, ifelse(rflw$variable  == "FLWtiller12", 2012, ifelse(rflw$variable  == "FLWtiller13", 2013, ifelse(rflw$variable  == "FLWtiller14", 2014, ifelse(rflw$variable  == "FLWtiller15", 2015, ifelse(rflw$variable  == "FLWtiller16", 2016, NA)))))))
# View(rflw)


rseedmerge_ss <- left_join(rspike, rseed, by = c( "plot", "pos", "tag", "Endo", 
                                              "Birth Year","year", "tillerid"))
# View(rseedmerge_ss)

rseedmerge_ssf <- merge(rseedmerge_ss, rflw, by = c("plot", "pos", "tag", "Endo", 
                                                    "Birth Year","year"), all = TRUE)
# View(rseedmerge_ssf)


# Pulling out the seed production estimates for the "Old" POAL data --------
pold_seed <- POAL_data_old %>%
  mutate("Birth Year" = year(`Date`)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate(seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2008", "seeds_InflA3", 'seeds_InflB3', 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

pold_seed$year<- ifelse(pold_seed$variable == "seed2008", 2008, ifelse(pold_seed$variable  == "seeds_InflA3", 2009, ifelse(pold_seed$variable == "seeds_InflB3", 2009, ifelse(pold_seed$variable == "seed2010", 2010, ifelse(pold_seed$variable  == "seed2011", 2011,ifelse(pold_seed$variable  == "seed2012", 2012, ifelse(pold_seed$variable  == "seed2013", 2013,ifelse(pold_seed$variable == "seed2014", 2014,ifelse(pold_seed$variable == "seed2015", 2015,ifelse(pold_seed$variable  == "seed2016", 2016, NA))))))))))


# View(pold_seed)

pold_spike <- POAL_data_old%>% 
  mutate("Birth Year" = year(`Date`)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("Spikelets_InflA1", "Spikelets_InflB1","spikelets_InflA3",
                       "spikelets_InflB3", "spikelets_InflA4","spikelets_InflB4",
                       "spikelets_InflC4","spikelets_InflA5", "spikelets_inflB5", 
                       "spikelets_inflC5", "spikelets_inflA6", 
                       "spikelets_inflA7", "spikelets_InflB7",  "spikelets_inflC7",
                       "spikelets_inflA8", "spikelets_InflB8", "spikelets_inflC8",
                       "spikelets_inflA9", "spikelets_InflB9","spikelets_inflC9",
                       "spikelets_inflA10", "spikelets_InflB10","spikelets_inflC10"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))
pold_spike$year<- ifelse(pold_spike$variable == "Spikelets_InflA1", 2008, ifelse(pold_spike$variable == "Spikelets_InflB1", 2008, ifelse(pold_spike$variable == "spikelets_InflA3", 2009, ifelse(pold_spike$variable  == "spikelets_InflB3", 2009, ifelse(pold_spike$variable  == "spikelets_InflA4", 2010, ifelse(pold_spike$variable  == "spikelets_InflB4", 2010, ifelse(pold_spike$variable  == "spikelets_InflC4", 2010, ifelse(pold_spike$variable  == "spikelets_InflA5", 2011, ifelse(pold_spike$variable  == "spikelets_inflB5", 2011, ifelse(pold_spike$variable  == "spikelets_inflC5", 2011,ifelse(pold_spike$variable == "spikelets_inflA6", 2012, ifelse(pold_spike$variable == "spikelets_inflA7", 2013, ifelse(pold_spike$variable == "spikelets_InflB7", 2013, ifelse(pold_spike$variable == "spikelets_inflC7", 2013, ifelse(pold_spike$variable  == "spikelets_inflA8", 2014, ifelse(pold_spike$variable == "spikelets_InflB8", 2014, ifelse(pold_spike$variable == "spikelets_inflC8", 2014, ifelse(pold_spike$variable == "spikelets_inflA9", 2015, ifelse(pold_spike$variable == "spikelets_InflB9", 2015, ifelse(pold_spike$variable ==   "spikelets_inflC9", 2015, ifelse(pold_spike$variable == "spikelets_inflA10", 2016, ifelse(pold_spike$variable == "spikelets_InflB10", 2016, ifelse(pold_spike$variable == "spikelets_inflC10", 2016, NA)))))))))))))))))))))))
# View(pold_spike)

# We already have the FlwTiller data within poldflw dataframe
poldflw <- POAL_data_old %>% 
  mutate("Birth Year" = year(`Date`)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("Flwtillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
poldflw$year<- ifelse(poldflw$variable == "Flwtillers1", 2008, ifelse(poldflw$variable  == "FlwTillers3", 2009, ifelse(poldflw$variable  == "FlwTillers4", 2010, ifelse(poldflw$variable  == "FlwTillers5", 2011, ifelse(poldflw$variable  == "FlwTillers6", 2012, ifelse(poldflw$variable  == "FlwTillers7", 2013,ifelse(poldflw$variable == "FlwTillers8", 2014,ifelse(poldflw$variable == "FlwTillers9", 2015,ifelse(poldflw$variable  == "FlwTillers10", 2016, NA)))))))))
# View(poldflw)

pold_seedmerge_ss <- left_join(pold_spike, pold_seed, by = c( "plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year", "tillerid"))
# View(pold_seedmerge_ss)

pold_seedmerge_ssf <- merge(pold_seedmerge_ss, poldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
# View(pold_seedmerge_ssf)



# Pulling out the seed, spikelet and flw tiller info from the POAL New Recruits
rold_seed <-POAL_data_old_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2009", "seed2010", "seed2011", 
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))
rold_seed$year<- ifelse(rold_seed$variable == "seed2009", 2009,ifelse(rold_seed$variable == "seed2010", 2010, ifelse(rold_seed$variable  == "seed2011", 2011, ifelse(rold_seed$variable  == "seed2012", 2012, ifelse(rold_seed$variable  == "seed2013", 2013, ifelse(rold_seed$variable  == "seed2014", 2014, ifelse(rold_seed$variable  == "seed2015", 2015, ifelse(rold_seed$variable  == "seed2016", 2016, NA))))))))

# View(rold_seed)

rold_spike <- POAL_data_old_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2009 = NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2010", "spike2011", "spikelets_inflA5", "spikelets_infl13","spike1","spike2", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("spike1",variable)~"A",
                              grepl("spike2", variable)~"B"))
rold_spike$year<- ifelse(rold_spike$variable == "spike2010", 2010, ifelse(rold_spike$variable == "spike2011", 2011, ifelse(rold_spike$variable == "spikelets_inflA5", 2012, ifelse(rold_spike$variable  == "spikelets_infl13", 2013, ifelse(rold_spike$variable  == "spike1", 2014, ifelse(rold_spike$variable == "spike2", 2014, ifelse(rold_spike$variable  == "spikelets_infl15", 2015, ifelse(rold_spike$variable  == "spikelets_infl16", 2016, NA))))))))
# View(rold_spike)

# We already have the FlwTiller data within roldflw dataframe
roldflw <- POAL_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWTiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
roldflw$year<- ifelse(roldflw$variable == "FLWTiller09", 2009, ifelse(roldflw$variable == "FLWtiller10", 2010, ifelse(roldflw$variable == "FLWtiller11", 2011, ifelse(roldflw$variable  == "FLWtiller12", 2012, ifelse(roldflw$variable  == "FLWtiller13", 2013, ifelse(roldflw$variable  == "FLWtiller14", 2014, ifelse(roldflw$variable  == "FLWtiller15", 2015, ifelse(roldflw$variable  == "FLWtiller16", 2016, NA))))))))
# View(roldflw)

rold_seedmerge_ss <- left_join(rold_spike, rold_seed, by = c( "plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year", "tillerid"))
# View(rold_seedmerge_ss)

rold_seedmerge_ssf <- merge(rold_seedmerge_ss, roldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
  

# View(rold_seedmerge_ssf)

pseedmerge_ssf <- pseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                  "tillerid", "flw","spikelets", "seed")]
pold_seedmerge_ssf <- pold_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                           "tillerid", "flw","spikelets", "seed")]
rseedmerge_ssf <- rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                   "tillerid", "flw","spikelets", "seed")]
rold_seedmerge_ssf <- rold_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                           "tillerid", "flw","spikelets", "seed")]


str(pseedmerge_ssf)
str(rseedmerge_ssf)
str(pold_seedmerge_ssf)
str(rold_seedmerge_ssf)

POALrepro <- pseedmerge_ssf %>% 
  rbind(rseedmerge_ssf) %>% 
  rbind(pold_seedmerge_ssf) %>% 
  rbind(rold_seedmerge_ssf) %>% 
  mutate(species = "POAL")
  


POALrepro <- POALrepro[!(is.na(POALrepro$flw)),]

# View(POALrepro)

# create version with t and t1 the will be used for merge with endo_demog_long
POALrepro_t1 <-POALrepro %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, flw_t1 = flw, seed_t1 = seed, spikelets_t1 = spikelets, tillerid_t1 = tillerid) %>%  
  mutate(year_t = year_t1 - 1)
# View(POAL_repro_t1)

POALrepro_t <-POALrepro %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, flw_t = flw, seed_t = seed, spikelets_t = spikelets, tillerid_t = tillerid) %>% 
  mutate(year_t1 = year_t + 1)
# View(POAL_repro_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
POALrepro_2 <- POALrepro_t1 %>% 
  full_join(POALrepro_t, by = c("plot", "pos", "tag", "Endo", 
                              "Birth Year", "year_t", "year_t1", "species"), all.x = all, all.y = all)
 # View(POALrepro_2)










# Combining seed productions measurements across years for the “New” POSY data -------------

## recoding for the year of measurement
## merging these measurements into one dataframe
po_seed <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2008", "seeds_InflA2", "seeds_InflB2", 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
         

po_seed$year<- ifelse(po_seed$variable == "seed2008", 2008,ifelse(po_seed$variable == "seeds_InflA2", 2009, ifelse(po_seed$variable  == "seeds_InflB2", 2009, ifelse(po_seed$variable  == "seed2010", 2010, ifelse(po_seed$variable  == "seed2011", 2011, ifelse(po_seed$variable  == "seed2012", 2012, ifelse(po_seed$variable  == "seed2013", 2013,ifelse(po_seed$variable == "seed2014", 2014,ifelse(po_seed$variable == "seed2015", 2015,ifelse(po_seed$variable  == "seed2016", 2016, NA))))))))))


# View(po_seed)

po_spike <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2008", "spikelets_InflA2", "spikelets_InflB2", "spikelets_inflA3", 
                       "spikelets_inflB3", "spikelets_inflA4","spikelets_inflB4",
                       "spikelets_inflA5", "spikelets_inflA6", "spikelets_InflB6", 
                       "spikelets_inflC6", "spikelets_inflA7", "spikelets_InflB7", 
                       "spikelets_inflC7","spikelets_inflA8", "spikelets_InflB8", 
                       "spikelets_inflC8"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))
po_spike$year<- ifelse(po_spike$variable == "spike2008", 2008,ifelse(po_spike$variable == "spikelets_InflA2", 2009, ifelse(po_spike$variable  == "spikelets_InflB2", 2009, ifelse(po_spike$variable  == "spikelets_inflA3", 2010, ifelse(po_spike$variable  == "spikelets_inflB3", 2010, ifelse(po_spike$variable  == "spikelets_inflB3", 2010, ifelse(po_spike$variable  == "spikelets_inflA4", 2011, ifelse(po_spike$variable == "spikelets_inflB4", 2011, ifelse(po_spike$variable == "spikelets_inflA5", 2012,ifelse(po_spike$variable == "spikelets_inflA6", 2013, ifelse(po_spike$variable == "spikelets_InflB6", 2013, ifelse(po_spike$variable == "spikelets_inflC6", 2013, ifelse(po_spike$variable  == "spikelets_inflA7", 2014, ifelse(po_spike$variable == "spikelets_InflB7", 2014, ifelse(po_spike$variable == "spikelets_inflC7", 2014, ifelse(po_spike$variable == "spikelets_inflA8", 2015, ifelse(po_spike$variable == "spikelets_InflB8", 2015, ifelse(po_spike$variable ==   "spikelets_inflC8", 2015, ifelse(po_spike$variable == "spikelets_inflA9", 2016, ifelse(po_spike$variable == "spikelets_InflB9", 2016, ifelse(po_spike$variable == "spikelets_inflC9", 2016, NA)))))))))))))))))))))
# View(po_spike)



po_flw <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8"), 
       value.name = "flw") 
po_flw$year<- ifelse(po_flw$variable == "Flwtillers1", 2008, ifelse(po_flw$variable  == "FlwTillers2", 2009, ifelse(po_flw$variable  == "FlwTillers3", 2010, ifelse(po_flw$variable  == "FlwTillers4", 2011, ifelse(po_flw$variable  == "FlwTillers5", 2012, ifelse(po_flw$variable  == "FlwTillers6", 2013,ifelse(po_flw$variable == "FlwTillers7", 2014,ifelse(po_flw$variable == "FlwTillers8", 2015,ifelse(po_flw$variable  == "FlwTillers9", 2016, NA)))))))))
# View(po_flw)

po_seedmerge_ss <- left_join(po_spike, po_seed, by = c("plot", "pos", "tag", "Endo", 
                                                   "Birth Year","year","tillerid"))
# View(po_seedmerge_ss)


po_seedmerge_ssf <- merge(po_seedmerge_ss,po_flw, by = c("plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year"))

# View(po_seedmerge_ssf)




# Combining repro measurements across years for the "New" POSY recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
po_rseed <- POSY_data_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2009", "seed2010", "seed2011",
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))


po_rseed$year<- ifelse(po_rseed$variable == "seed", 2009, ifelse(po_rseed$variable  == "seed2010", 2010, ifelse(po_rseed$variable  == "seed2011", 2011, ifelse(po_rseed$variable  == "seed2012", 2012, ifelse(po_rseed$variable  == "seed2013", 2013,ifelse(po_rseed$variable == "seed2014", 2014,ifelse(po_rseed$variable == "seed2015", 2015,ifelse(po_rseed$variable  == "seed2016", 2016, NA))))))))


# View(po_rseed)

po_rspike <- POSY_data_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike_2009 = NA, spike_2010 = NA, spike_2011 = NA,) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike_2009", "spike_2010", "spike_2011", 
                       "spikelets_inflA5", "spikelets_inflA5__1", 
                       "spikelets_infl14","spike1", "spike2", 
                       "spike3", "spike1__1", "spike2__1", "spike3__1"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                              grepl("spike2", variable) ~ "B",
                              grepl("spike3", variable) ~ "C"))
po_rspike$year<- ifelse(po_rspike$variable  == "spike_2009", 2009, ifelse(po_rspike$variable  == "spike_2010", 2010, ifelse(po_rspike$variable == "spike_2011", 2011, ifelse(po_rspike$variable  == "spikelets_inflA5", 2012, ifelse(po_rspike$variable == "spikelets_inflA5__1", 2013, ifelse(po_rspike$variable == "spikelets_InflB6", 2013, ifelse(po_rspike$variable == "spikelets_inflC6", 2013, ifelse(po_rspike$variable  == "spikelets_infl14", 2014, ifelse(po_rspike$variable == "spike1", 2015, ifelse(po_rspike$variable == "spike2", 2015, ifelse(po_rspike$variable ==   "spike3", 2015, ifelse(po_rspike$variable == "spike1__1", 2016, ifelse(po_rspike$variable == "spike2__1", 2016, ifelse(po_rspike$variable == "spike3__1", 2016, NA))))))))))))))
# View(po_rspike)



po_rflw <- POSY_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(FLWtiller09 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c( "FLWtiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_rflw$year<- ifelse(po_rflw$variable == "FLWtiller09", 2009, ifelse(po_rflw$variable == "FLWtiller10", 2010, ifelse(po_rflw$variable == "FLWtiller11", 2011, ifelse(po_rflw$variable  == "FLWtiller12", 2012, ifelse(po_rflw$variable  == "FLWtiller13", 2013, ifelse(po_rflw$variable  == "FLWtiller14", 2014, ifelse(po_rflw$variable  == "FLWtiller15", 2015, ifelse(po_rflw$variable  == "FLWtiller16", 2016, NA))))))))
po_rflw$flw <- ifelse(grepl("09", po_rflw$`Birth Year`) & grepl("09", po_rflw$variable), 0, NA)
# View(po_rflw)


po_rmerge_ss <- left_join(po_rspike, po_rseed, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_rmerge_ss)

po_rmerge_ssf <- merge(po_rmerge_ss, po_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_rmerge_ssf)






# Combining repro measurements across years for the “Old” POSY data ------------------


## Combining data across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
po_oldsurv <- POSY_data_old %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"),
       measure.var = c("survive1", "survive3", "survive4", 
                       "survive5", "Survive6", "survive7", 
                       "survive8", "survive9", "survive10"),
       value.name = "surv") 
po_oldsurv$year<- ifelse(po_oldsurv$variable == "survive1", 2008, ifelse(po_oldsurv$variable  == "survive3", 2009, ifelse(po_oldsurv$variable  == "survive4", 2010, ifelse(po_oldsurv$variable  == "survive5", 2011, ifelse(po_oldsurv$variable  == "Survive6", 2012, ifelse(po_oldsurv$variable  == "survive7", 2013,ifelse(po_oldsurv$variable == "survive8", 2014,ifelse(po_oldsurv$variable == "survive9", 2015,ifelse(po_oldsurv$variable  == "survive10", 2016, NA)))))))))
# View(po_oldsurv)

po_oldgrow <- POSY_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("TotTillers1", "TotTillers3","TotTillers4", 
                       "TotTillers5", "TotTillers6","TotTillers7", 
                       "TotTillers8", "TotTillers9", "TotTillers10"), 
       value.name = "size") 
po_oldgrow$year<- ifelse(po_oldgrow$variable == "TotTillers1", 2008, ifelse(po_oldgrow$variable  == "TotTillers3", 2009, ifelse(po_oldgrow$variable  == "TotTillers4", 2010, ifelse(po_oldgrow$variable  == "TotTillers5", 2011, ifelse(po_oldgrow$variable  == "TotTillers6", 2012, ifelse(po_oldgrow$variable  == "TotTillers7", 2013, ifelse(po_oldgrow$variable == "TotTillers8", 2014, ifelse(po_oldgrow$variable == "TotTillers9", 2015, ifelse(po_oldgrow$variable  == "TotTillers10", 2016, NA)))))))))
# View(po_oldgrow)

po_oldflw <- POSY_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n",
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("FlwTillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
po_oldflw$year<- ifelse(po_oldflw$variable == "FlwTillers1", 2008, ifelse(po_oldflw$variable  == "FlwTillers3", 2009, ifelse(po_oldflw$variable  == "FlwTillers4", 2010, ifelse(po_oldflw$variable  == "FlwTillers5", 2011, ifelse(po_oldflw$variable  == "FlwTillers6", 2012, ifelse(po_oldflw$variable  == "FlwTillers7", 2013,ifelse(po_oldflw$variable == "FlwTillers8", 2014,ifelse(po_oldflw$variable == "FlwTillers9", 2015,ifelse(po_oldflw$variable  == "FlwTillers10", 2016, NA)))))))))
# View(po_oldflw)


po_oldmerge_sg <- merge(po_oldsurv, po_oldgrow, by = c( "plot","pos", "tag", "Endo", 
                                                        "Loc'n", "Birth Year", "TRT",
                                                        "Plant", "year"))
# View(po_oldmerge_sg)

po_oldmerge_sgf <- merge(po_oldmerge_sg, po_oldflw, by = c( "plot","pos", "tag", "Endo", 
                                                            "Loc'n", "Birth Year", "TRT",
                                                            "Plant", "year"))
# View(po_ldmerge_sgf)

# getting a dataframe with t and t_1
po_oldmerge_t1 <-po_oldmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(po_oldmerge_t1)

po_oldmerge_t <-po_oldmerge_sgf %>%
  filter(year != max(year)) %>% 
  select(-surv) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(po_oldmerge_t)

po_oldmerge <- po_oldmerge_t1 %>% 
  full_join(po_oldmerge_t, by = c("plot","pos", "tag", "Endo", 
                                  "Loc'n", "Birth Year", "TRT",
                                  "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(po_oldmerge)


# Combining repro measurements across years for the “Old” POSY recruits data --------

## Combining data for recruits across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
po_roldsurv <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive09 = "NA") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("survive09", "survive10", "survive11", "survive12", "survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
po_roldsurv$year<- ifelse(po_roldsurv$variable == "survive09", 2009, ifelse(po_roldsurv$variable == "survive10", 2010, ifelse(po_roldsurv$variable == "survive11", 2011, ifelse(po_roldsurv$variable  == "survive12", 2012, ifelse(po_roldsurv$variable  == "survive13", 2013, ifelse(po_roldsurv$variable  == "Survive14", 2014, ifelse(po_roldsurv$variable  == "Survive15", 2015, ifelse(po_roldsurv$variable  == "Survive16", 2016, NA))))))))
# View(po_roldsurv)


po_roldgrow <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller09", "TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
po_roldgrow$year<- ifelse(po_roldgrow$variable == "TOTtiller09", 2009, ifelse(po_roldgrow$variable == "TOTtiller10", 2010, ifelse(po_roldgrow$variable == "TOTtiller11", 2011, ifelse(po_roldgrow$variable  == "TOTtiller12", 2012, ifelse(po_roldgrow$variable  == "TOTtiller13", 2013, ifelse(po_roldgrow$variable  == "TOTtiller14", 2014, ifelse(po_roldgrow$variable  == "TOTtiller15", 2015, ifelse(po_roldgrow$variable  == "TOTtiller16", 2016, NA))))))))
# View(po_roldgrow)

po_roldflw <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(FLWtiller09 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller09", "FLWtiller10","FLWtiller11", "FLWTiller12", 
                       "FLWTiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_roldflw$year<- ifelse(po_roldflw$variable == "FLWtiller09", 2009, ifelse(po_roldflw$variable == "FLWtiller10", 2010, ifelse(po_roldflw$variable == "FLWtiller11", 2011, ifelse(po_roldflw$variable  == "FLWTiller12", 2012, ifelse(po_roldflw$variable  == "FLWTiller13", 2013, ifelse(po_roldflw$variable  == "FLWtiller14", 2014, ifelse(po_roldflw$variable  == "FLWtiller15", 2015, ifelse(po_roldflw$variable  == "FLWtiller16", 2016, NA))))))))
# View(po_roldflw)

po_roldseedmerge_ss <- left_join(po_roldspike, po_roldseed, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year", "tillerid"))
# View(po_roldmerge_ss)

po_roldseedmerge_ssf <- merge(po_roldseedmerge_ss, po_roldflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"), all = TRUE)
# View(po_roldmerge_ssf)






# Combining the old and new and original and recruit POSY repro dataframes ---------
po_seedmerge_ssf <- po_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                   "tillerid", "flw","spikelets", "seed")]
po_oldseedmerge_ssf <- po_oldseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                           "tillerid", "flw","spikelets", "seed")]
po_rseedmerge_ssf <- po_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                   "tillerid", "flw","spikelets", "seed")]
po_roldseedmerge_ssf <- po_roldseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                           "tillerid", "flw","spikelets", "seed")]


str(po_seedmerge_ssf)
str(po_rseedmerge_ssf)
str(po_oldseedmerge_ssf)
str(po_roldseedmerge_ssf)

POSYrepro <- po_seedmerge_ssf %>% 
  rbind(po_rseedmerge_ssf) %>% 
  rbind(po_oldseedmerge_ssf) %>% 
  rbind(po_roldseedmerge_ssf) %>% 
  mutate(species = "POSY")

POSYrepro <- POSYrepro[!(is.na(POSYrepro$flw)),]

# View(POSYrepro)

# create version with t and t1 the will be used for merge with endo_demog_long
POSYrepro_t1 <-POSYrepro %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, flw_t1 = flw, seed_t1 = seed, spikelets_t1 = spikelets, tillerid_t1 = tillerid) %>%  
  mutate(year_t = year_t1 - 1)
# View(POSYrepro_t1)

POSYrepro_t <-POSYrepro %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, flw_t = flw, seed_t = seed, spikelets_t = spikelets, tillerid_t = tillerid) %>% 
  mutate(year_t1 = year_t + 1)
# View(POSYrepro_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
POSYrepro_2 <- POSYrepro_t1 %>% 
  merge(POSYrepro_t, by = c("plot", "pos", "tag", "Endo", 
                                "Birth Year", "year_t", "year_t1", "species"))
# View(POSYrepro_2)












# Combining repro measurements across years for the LOAR data ------------------

## Combining measurements across years for reproduction measurements

lseed <- LOAR_data %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG" ) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed1", "seed2", "seed3","seed5",
                       "seed6(2012)", "seed7(2013)", "seed8(2014)", 
                       "seed9(2015)", "seed10(2016)"),
       value.name = "seed") 
lseed$year<- ifelse(lseed$variable == "seed1", 2008, ifelse(lseed$variable  == "seed2", 2009, ifelse(lseed$variable  == "seed3", 2010, ifelse(lseed$variable  == "seed5", 2011, ifelse(lseed$variable  == "seed6(2012)", 2012, ifelse(lseed$variable  == "seed7(2013)", 2013,ifelse(lseed$variable == "seed8(2014)", 2014,ifelse(lseed$variable == "seed9(2015)", 2015,ifelse(lseed$variable  == "seed10(2016)", 2016, NA)))))))))
View(lseed)

lseed2008 <- LOAR_data_seed2008 %>% 
  rename("tag" = "Tag" ) %>%
  distinct(tag) 
  melt(id.var = c("tag"),
       measure.var = c("seed1", "seed2", "seed3","seed5",
                       "seed6(2012)", "seed7(2013)", "seed8(2014)", 
                       "seed9(2015)", "seed10(2016)"),
       value.name = "seed") 
lseed2008$year<- ifelse(lseed2008$variable == "seed1", 2008, ifelse(lseed2008$variable  == "seed2", 2009, ifelse(lseed2008$variable  == "seed3", 2010, ifelse(lseed2008$variable  == "seed5", 2011, ifelse(lseed2008$variable  == "seed6(2012)", 2012, ifelse(lseed2008$variable  == "seed7(2013)", 2013,ifelse(lseed2008$variable == "seed8(2014)", 2014,ifelse(lseed2008$variable == "seed9(2015)", 2015,ifelse(lseed2008$variable  == "seed10(2016)", 2016, NA)))))))))
View(lseed2008)

lseed2009 <- LOAR_data_seed2009%>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG" ) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed1", "seed2", "seed3","seed5",
                       "seed6(2012)", "seed7(2013)", "seed8(2014)", 
                       "seed9(2015)", "seed10(2016)"),
       value.name = "seed") 
lseed2009$year<- ifelse(lseed2009$variable == "seed1", 2008, ifelse(lseed2009$variable  == "seed2", 2009, ifelse(lseed2009$variable  == "seed3", 2010, ifelse(lseed2009$variable  == "seed5", 2011, ifelse(lseed2009$variable  == "seed6(2012)", 2012, ifelse(lseed2009$variable  == "seed7(2013)", 2013,ifelse(lseed2009$variable == "seed8(2014)", 2014,ifelse(lseed2009$variable == "seed9(2015)", 2015,ifelse(lseed2009$variable  == "seed10(2016)", 2016, NA)))))))))
View(lseed2009)

lspike <- LOAR_data %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("TOTtillers1", "TOTtiller2", "TOTtillers3",
                       "TOTtillers5", "TOTtillers6","TOTtillers7", 
                       "TOTtillers8", "TOTtillers9", "TOTtillers10"), 
       value.name = "size") 
lspike$year<- ifelse(lspike$variable == "TOTtillers1", 2008, ifelse(lspike$variable  == "TOTtiller2", 2009, ifelse(lspike$variable  == "TOTtillers3", 2010, ifelse(lspike$variable  == "TOTtillers5", 2011, ifelse(lspike$variable  == "TOTtillers6", 2012, ifelse(lspike$variable  == "TOTtillers7", 2013, ifelse(lspike$variable == "TOTtillers8", 2014, ifelse(lspike$variable == "TOTtillers9", 2015, ifelse(lspike$variable  == "TOTtillers10", 2016, NA)))))))))
# View(lspike)

lflw <- LOAR_data %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("FlwTillers1", "FLWTiller2", "FLWTillers3", 
                       "FLWTillers5", "FLWTillers6", "FLWTillers7", 
                       "FLWTillers8", "FLWTillers9", "FLWTillers10"), 
       value.name = "flw") 
lflw$year<- ifelse(lflw$variable == "FlwTillers1", 2008, ifelse(lflw$variable  == "FLWTiller2", 2009, ifelse(lflw$variable  == "FLWTillers3", 2010, ifelse(lflw$variable  == "FLWTillers5", 2011, ifelse(lflw$variable  == "FLWTillers6", 2012, ifelse(lflw$variable  == "FLWTillers7", 2013,ifelse(lflw$variable == "FLWTillers8", 2014,ifelse(lflw$variable == "FLWTillers9", 2015,ifelse(lflw$variable  == "FLWTillers10", 2016, NA)))))))))
# View(lflw)

l_seedmerge_ss <- merge(lseed, lspike, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year", "tillerid"))
# View(l_merge_sg)

l_seedmerge_ssf <- merge(l_seedmerge_ss, lflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(l_merge_sgf)

# getting a dataframe with t and t_1
l_merge_t1 <-l_merge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, seed_t1 = seed, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(l_merge_t1)

l_merge_t <-l_merge_sgf %>%
  filter(year != max(year)) %>% 
  select(-seed) %>% 
  rename(year_t = year, size_t = size, flw_t = flw) 
# View(l_merge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
l_merge <- l_merge_t1 %>% 
  full_join(l_merge_t, by = c("plot", "pos", "tag", "Endo", 
                              "Loc'n", "Birth Year", "TRT",
                              "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(l_merge)


# Combining repro measurements across years for the LOAR recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
l_rseed <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  mutate(seedive10 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("seedive10", "seedive11", "seedive12", "seedive13", "seedive14", 
                       "seedive15", "seedive16"),
       value.name = "seed") 
l_rseed$year<- ifelse(l_rseed$variable == "seedive10", 2010, ifelse(l_rseed$variable == "seedive11", 2011, ifelse(l_rseed$variable  == "seedive12", 2012, ifelse(l_rseed$variable  == "seedive13", 2013, ifelse(l_rseed$variable  == "seedive14", 2014, ifelse(l_rseed$variable  == "seedive15", 2015, ifelse(l_rseed$variable  == "seedive16", 2016, NA)))))))
# View(l_rseed)

l_rspike <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TotTillers10","Tot11", "Tot12", 
                       "Tot13", "Tot14", "Tot15", 
                       "Tot16"),
       value.name = "size") 
l_rspike$year<- ifelse(l_rspike$variable == "TotTillers10", 2010, ifelse(l_rspike$variable == "Tot11", 2011, ifelse(l_rspike$variable  == "Tot12", 2012, ifelse(l_rspike$variable  == "Tot13", 2013, ifelse(l_rspike$variable  == "Tot14", 2014, ifelse(l_rspike$variable  == "Tot15", 2015, ifelse(l_rspike$variable  == "Tot16", 2016, NA)))))))
# View(l_rspike)

l_rflw <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FlwTillers10","Flw11", "Flw12", 
                       "Flw13", "Flw14", "Flw15",
                       "Flw16"),
       value.name = "flw") 
l_rflw$year<- ifelse(l_rflw$variable == "FlwTillers10", 2010, ifelse(l_rflw$variable == "Flw11", 2011, ifelse(l_rflw$variable  == "Flw12", 2012, ifelse(l_rflw$variable  == "Flw13", 2013, ifelse(l_rflw$variable  == "Flw14", 2014, ifelse(l_rflw$variable  == "Flw15", 2015, ifelse(l_rflw$variable  == "Flw16", 2016, NA)))))))
# View(l_rflw)

l_rmerge_sg <- merge(l_rseed, l_rspike, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(l_rmerge_sg)
=======
po_roldmerge_sg <- merge(po_roldsurv, po_roldgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_roldmerge_sg)
>>>>>>> 3b195afb74218a5ee8d92add3f8b76cff775ad42

po_roldmerge_sgf <- merge(po_roldmerge_sg, po_roldflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(po_roldmerge_sgf)

## getting a dataframe with time t and t_1
po_roldmerge_t1 <-po_roldmerge_sgf %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(po_roldmerge_t1)

po_roldmerge_t <-po_roldmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(po_roldmerge_t)

po_roldmerge <- po_roldmerge_t1 %>% 
  full_join(po_roldmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(po_roldmerge)







# Combining the old and new and original and recruit POSY dataframes ---------
po_merge <- po_merge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                       "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                       "year_t", "size_t", "flw_t")]
po_ldmerge <- po_oldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                            "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                            "year_t", "size_t", "flw_t")]
po_rmerge <- po_rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "size_t", "flw_t")]
po_roldmerge <- po_roldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                               "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                               "year_t", "size_t", "flw_t")]

POSY <- po_merge %>% 
  rbind(po_oldmerge) %>% 
  rbind(po_roldmerge) %>% 
  rbind(po_rmerge) %>% 
  mutate(species = "POSY") %>% 
  mutate(surv_t1 = ifelse(surv_t1 == 'NA', NA, surv_t1))
POSY <- POSY[!(is.na(POSY$surv_t1)),]

View(POSY)









# What we want our final data frame to look like
LTREB_endodemog <- data.frame(colnames("quad", "species", "origin", "plot", "pos", "id", "surv_t1", "size_t1", "seed_t1", "flower_t1", "size_t", "endo", "birth", "year_t", "year_t1"))






=======
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
LTREB_endodemog <- 
  read.csv("~/Documents/R projects/LTREBendodemog/endo_demog_long.csv")

View(LTREB_endodemog)
str(LTREB_endodemog)
dim(LTREB_endodemog)
>>>>>>> 92e3a4d30db8bc4df7a7d6315422ba4a8dd9c6f9

