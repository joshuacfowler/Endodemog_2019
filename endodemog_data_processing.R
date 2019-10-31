## Authors: Josh and Tom	## Grass endophyte population model
## Purpose: Create a script that imports Endodemog data, perform all raw data manipulation to set up data lists for Survival and Growth models and for the flowering tiller and seed production models,	
## and create an .RData object that can be loaded for analysis	
## Last Update: Jul 20, 2019
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)


##############################################################################
####### Reading in endo_demog_long ------------------------------
##############################################################################


# This is the main compiled data sheet that we will merge with the flowering and seed data
LTREB_endodemog <- 
  read.csv(file = "~/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t = na_if(size_t, 0)) %>% 
  mutate(size_t1 = na_if(size_t1, 0)) %>%  
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>%
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
  mutate(species = case_when(species == "ELVI" & plot == 101 ~ "ELRI", species == "ELVI" & plot != 101 ~ "ELVI",  # This is for the compiled data where a ELRI data point in plot 101, tag 2004 is labelled as ELVI
                                   species == "ELRI" ~ "ELRI",
                                   species == "FESU" ~ "FESU",
                                   species == "AGPE" ~ "AGPE",
                                   species == "POAL" ~ "POAL",
                                   species == "POSY" ~ "POSY",
                                   species == "LOAR" ~ "LOAR"))  %>%    
  mutate(species_index = as.integer(recode_factor(species,                   
                                                  "AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                                                  "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                                                  "POSY" = 7))) %>% 
  mutate(year_t_index = as.integer(recode(year_t, 
                                          '2007' = 1, '2008' = 2, '2009' = 3, 
                                          '2010' = 4, '2011' = 5, '2012' = 6, 
                                          '2013' = 7, '2014' = 8, '2015' = 9, 
                                          '2016' = 10, '2017' = 11))) %>%             
  mutate(year_t1_index = as.integer(recode(year_t1, 
                                           '2008' = 2, '2009' = 3, '2010' = 4, 
                                           '2011' = 5, '2012' = 6, '2013' = 7, 
                                           '2014' = 8, '2015' = 9, '2016' = 10, 
                                           '2017' = 11, '2018' = 12))) %>%               
  mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, 
                                          origin == "R" ~ 1, 
                                          origin != "R" | origin != "O" ~ 1))) %>%
  mutate(plot_fixed = as.integer(case_when(species == "LOAR" & id == "39_1B" ~ "39", # This is for a copy error with LOAR individual 39_1B, which assigned it to plots 41-44
                                 plot != "R" ~ as.character(plot), 
                                 plot == "R" ~ as.character(origin)))) %>% 
  mutate(surv_t1 = as.integer(case_when(surv_t1 == 1 ~ 1,
                                   surv_t1 == 0 ~ 0,
                                   is.na(surv_t1) & birth == year_t1 ~ 1))) %>% 
  filter(duplicated(.) == FALSE)
# dim(LTREB_data)



########################################################################################################################
###### Reading in raw excel files which have the detailed Reproduction data.---------------------------
########################################################################################################################
# read in raw data from POAL
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
LOAR_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR")
LOAR_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR recruits")
LOAR_data_seed2008 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR seeds 2008", skip=1)
LOAR_data_seed2009 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR seeds 2009")

# Read in data from FESU
FESU_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU")
FESU_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU recruits")

# Read in data from ELVI
ELVI_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI(IN) originals up to 2016.xlsx", sheet = "ELVI")
ELVI_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVI recruits")
ELVI_data_seed <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVISeeds")
ELVI_data_r <- ELVI_data_r %>% 
  mutate(tag = paste(Plot, RecruitNo, sep = "_")) 
  

# Read in data from ELRI
ELRI_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI")
ELRI_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI recruits")
ELRI_data_seed2009 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI infl 09")
ELRI_data_r <- ELRI_data_r %>% 
  mutate(tag = paste(PLOT, RecruitNo, sep = "_")) 
# Read in data from AGPE
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



#################################################################################################################################################################
# Pulling out the seed production estimates. These are not measured for all plants, and so will go into a separate dataframe------------------------------
#################################################################################################################################################################
# Pulling out the seed production estimates for the "New" POAL data --------
pseed <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA2", "seeds_InflB2", 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))

pseed$year<- ifelse(pseed$variable == "seed2007", 2007, ifelse(pseed$variable == "seed2008", 2008,ifelse(pseed$variable == "seeds_InflA2", 2009, ifelse(pseed$variable  == "seeds_InflB2", 2009, ifelse(pseed$variable  == "seed2010", 2010, ifelse(pseed$variable  == "seed2011", 2011, ifelse(pseed$variable  == "seed2012", 2012, ifelse(pseed$variable  == "seed2013", 2013,ifelse(pseed$variable == "seed2014", 2014,ifelse(pseed$variable == "seed2015", 2015,ifelse(pseed$variable  == "seed2016", 2016, NA)))))))))))

pseed1 <- pseed %>% 
  filter(!is.na(seed))
 # View(pseed1)

pspike <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(spike2007 = NA, spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spikelets_InflA2", "spikelets_InflB2", "spikelets_inflA3", 
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
pspike$year<- ifelse(pspike$variable == "spike2007", 2007, ifelse(pspike$variable == "spike2008", 2008,ifelse(pspike$variable == "spikelets_InflA2", 2009, ifelse(pspike$variable  == "spikelets_InflB2", 2009, ifelse(pspike$variable  == "spikelets_inflA3", 2010, ifelse(pspike$variable  == "spikelets_inflB3", 2010, ifelse(pspike$variable  == "spikelets_inflC3", 2010, ifelse(pspike$variable  == "spikelets_inflA4", 2011,ifelse(pspike$variable == "spikelets_inflA5", 2012,ifelse(pspike$variable == "spikelets_inflA6", 2013, ifelse(pspike$variable == "spikelets_InflB6", 2013, ifelse(pspike$variable == "spikelets_inflC6", 2013, ifelse(pspike$variable  == "spikelets_inflA7", 2014, ifelse(pspike$variable == "spikelets_InflB7", 2014, ifelse(pspike$variable == "spikelets_inflC7", 2014, ifelse(pspike$variable == "spikelets_inflA8", 2015, ifelse(pspike$variable == "spikelets_InflB8", 2015, ifelse(pspike$variable ==   "spikelets_inflC8", 2015, ifelse(pspike$variable == "spikelets_inflA9", 2016, ifelse(pspike$variable == "spikelets_InflB9", 2016, ifelse(pspike$variable == "spikelets_inflC9", 2016, NA)))))))))))))))))))))
# View(pspike)
pspike1 <- pspike %>%  
  filter(!is.na(spikelets))
# View(pseed1)

pflw <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
pflw$year<- ifelse(pflw$variable == "flw2007", 2007, ifelse(pflw$variable == "Flwtillers1", 2008, ifelse(pflw$variable  == "FlwTillers2", 2009, ifelse(pflw$variable  == "FlwTillers3", 2010, ifelse(pflw$variable  == "FlwTillers4", 2011, ifelse(pflw$variable  == "FlwTillers5", 2012, ifelse(pflw$variable  == "FlwTillers6", 2013,ifelse(pflw$variable == "FlwTillers7", 2014,ifelse(pflw$variable == "FlwTillers8", 2015,ifelse(pflw$variable  == "FlwTillers9", 2016, NA))))))))))

# View(pflw)

pseedmerge_ss <- full_join(pspike1, pseed1, by = c( "plot", "pos", "tag", "Endo", 
                                         "Birth Year","year", "tillerid"))
# View(pseedmerge_ss)

pseedmerge_ssf <- merge(pseedmerge_ss, pflw, by = c("plot", "pos", "tag", "Endo", 
                                                    "Birth Year","year"), all = TRUE)
# View(pseedmerge_ssf)



# Pulling out the seed, spikelet and flw tiller info from the POAL New Recruits
rseed <-POAL_data_r %>% 
  mutate("Loc'n" = NA, "TRT" = NA) %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010","seed2011", "seed2012", "seed2013", "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))
rseed$year<- ifelse(rseed$variable == "seed2007", 2007, ifelse(rseed$variable == "seed2008", 2008, ifelse(rseed$variable == "seed2009", 2009, ifelse(rseed$variable == "seed2010", 2010, ifelse(rseed$variable  == "seed2011", 2011, ifelse(rseed$variable  == "seed2012", 2012, ifelse(rseed$variable  == "seed2013", 2013, ifelse(rseed$variable  == "seed2014", 2014, ifelse(rseed$variable  == "seed2015", 2015, ifelse(rseed$variable  == "seed2016", 2016, NA))))))))))
rseed1 <- rseed %>% 
  filter(!is.na(seed))
# View(rseed1)

rspike <- POAL_data_r %>% 
  mutate("Loc'n" = NA, "TRT" = NA) %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", "spikelets_inflA5", "spikelets_inflA5__1", "spikelets_infl14", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",))
rspike$year<- ifelse(rspike$variable == "spike2007", 2007, ifelse(rspike$variable == "spike2008", 2008, ifelse(rspike$variable == "spike2009", 2009, ifelse(rspike$variable == "spike2010", 2010, ifelse(rspike$variable == "spike2011", 2011, ifelse(rspike$variable == "spikelets_inflA5", 2012, ifelse(rspike$variable  == "spikelets_inflA5__1", 2013, ifelse(rspike$variable  == "spikelets_infl14", 2014, ifelse(rspike$variable  == "spikelets_infl15", 2015, ifelse(rspike$variable  == "spikelets_infl16", 2016, NA))))))))))
rspike1 <- rspike %>% 
  filter(!is.na(spikelets))
# View(rspike1)

# We already have the FlwTiller data within rflw dataframe
rflw <- POAL_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
rflw$year<- ifelse(rflw$variable == "flw2007", 2007, ifelse(rflw$variable == "flw2008", 2008, ifelse(rflw$variable == "flw2009", 2009, ifelse(rflw$variable == "FLWtiller10", 2010, ifelse(rflw$variable == "FLWtiller11", 2011, ifelse(rflw$variable  == "FLWtiller12", 2012, ifelse(rflw$variable  == "FLWtiller13", 2013, ifelse(rflw$variable  == "FLWtiller14", 2014, ifelse(rflw$variable  == "FLWtiller15", 2015, ifelse(rflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(rflw)


rseedmerge_ss <- full_join(rspike1, rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                              "Birth Year","year", "tillerid"))
# View(rseedmerge_ss)

rseedmerge_ssf <- merge(rseedmerge_ss, rflw, by = c("plot", "pos", "tag", "Endo", 
                                                    "Birth Year","year"), all = TRUE)
# View(rseedmerge_ssf)


# Pulling out the seed production estimates for the "Old" POAL data --------
pold_seed <- POAL_data_old %>%
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA3", 'seeds_InflB3', 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

pold_seed$year<- ifelse(pold_seed$variable == "seed2007", 2007, ifelse(pold_seed$variable == "seed2008", 2008, ifelse(pold_seed$variable  == "seeds_InflA3", 2009, ifelse(pold_seed$variable == "seeds_InflB3", 2009, ifelse(pold_seed$variable == "seed2010", 2010, ifelse(pold_seed$variable  == "seed2011", 2011,ifelse(pold_seed$variable  == "seed2012", 2012, ifelse(pold_seed$variable  == "seed2013", 2013,ifelse(pold_seed$variable == "seed2014", 2014,ifelse(pold_seed$variable == "seed2015", 2015,ifelse(pold_seed$variable  == "seed2016", 2016, NA)))))))))))

pold_seed1 <- pold_seed %>% 
  filter(!is.na(seed))

# View(pold_seed1)

pold_spike <- POAL_data_old%>% 
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate(spike2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "Spikelets_InflA1", "Spikelets_InflB1","spikelets_InflA3",
                       "spikelets_InflB3", "spikelets_InflA4","spikelets_InflB4",
                       "spikelets_InflC4","spikelets_InflA5", "spikelets_inflB5", 
                       "spikelets_inflC5", "spikelets_inflA6", 
                       "spikelets_inflA7", "spikelets_InflB7",  "spikelets_inflC7",
                       "spikelets_inflA8", "spikelets_InflB8", "spikelets_inflC8",
                       "spikelets_inflA9", "spikelets_InflB9","spikelets_inflC9",
                       "spikelets_inflA10", "spikelets_InflB10","spikelets_inflC10"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("InflA", variable) ~ "A", 
                              grepl("InflB", variable) ~ "B",
                              grepl("InflC", variable) ~ "C",
                              grepl("inflA", variable) ~ "A", 
                              grepl("inflB", variable) ~ "B",
                              grepl("inflC", variable) ~ "C"))
pold_spike$year<- ifelse(pold_spike$variable == "spike2007", 2007, ifelse(pold_spike$variable == "Spikelets_InflA1", 2008, ifelse(pold_spike$variable == "Spikelets_InflB1", 2008, ifelse(pold_spike$variable == "spikelets_InflA3", 2009, ifelse(pold_spike$variable  == "spikelets_InflB3", 2009, ifelse(pold_spike$variable  == "spikelets_InflA4", 2010, ifelse(pold_spike$variable  == "spikelets_InflB4", 2010, ifelse(pold_spike$variable  == "spikelets_InflC4", 2010, ifelse(pold_spike$variable  == "spikelets_InflA5", 2011, ifelse(pold_spike$variable  == "spikelets_inflB5", 2011, ifelse(pold_spike$variable  == "spikelets_inflC5", 2011,ifelse(pold_spike$variable == "spikelets_inflA6", 2012, ifelse(pold_spike$variable == "spikelets_inflA7", 2013, ifelse(pold_spike$variable == "spikelets_InflB7", 2013, ifelse(pold_spike$variable == "spikelets_inflC7", 2013, ifelse(pold_spike$variable  == "spikelets_inflA8", 2014, ifelse(pold_spike$variable == "spikelets_InflB8", 2014, ifelse(pold_spike$variable == "spikelets_inflC8", 2014, ifelse(pold_spike$variable == "spikelets_inflA9", 2015, ifelse(pold_spike$variable == "spikelets_InflB9", 2015, ifelse(pold_spike$variable ==   "spikelets_inflC9", 2015, ifelse(pold_spike$variable == "spikelets_inflA10", 2016, ifelse(pold_spike$variable == "spikelets_InflB10", 2016, ifelse(pold_spike$variable == "spikelets_inflC10", 2016, NA))))))))))))))))))))))))

pold_spike1 <- pold_spike %>% 
  filter(!is.na(spikelets))
# View(pold_spike1)

poldflw <- POAL_data_old %>% 
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
poldflw$year<- ifelse(poldflw$variable == "flw2007", 2007, ifelse(poldflw$variable == "Flwtillers1", 2008, ifelse(poldflw$variable  == "FlwTillers3", 2009, ifelse(poldflw$variable  == "FlwTillers4", 2010, ifelse(poldflw$variable  == "FlwTillers5", 2011, ifelse(poldflw$variable  == "FlwTillers6", 2012, ifelse(poldflw$variable  == "FlwTillers7", 2013,ifelse(poldflw$variable == "FlwTillers8", 2014,ifelse(poldflw$variable == "FlwTillers9", 2015,ifelse(poldflw$variable  == "FlwTillers10", 2016, NA))))))))))
# View(poldflw)

pold_seedmerge_ss <- full_join(pold_spike1, pold_seed1, by = c( "plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year", "tillerid"))
# View(pold_seedmerge_ss)

pold_seedmerge_ssf <- merge(pold_seedmerge_ss, poldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
# View(pold_seedmerge_ssf)



# Pulling out the seed, spikelet and flw tiller info from the POAL New Recruits
rold_seed <-POAL_data_old_r %>% 
  rename("tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate("Birth Year" = case_when(Date == "Q2010" ~ "2010",
                                  Date != "Q2010" ~ Date)) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", 
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B")) 
rold_seed$year<- ifelse(rold_seed$variable == "seed2007", 2007, ifelse(rold_seed$variable == "seed2008", 2008, ifelse(rold_seed$variable == "seed2008", 2008, ifelse(rold_seed$variable == "seed2009", 2009,ifelse(rold_seed$variable == "seed2010", 2010, ifelse(rold_seed$variable  == "seed2011", 2011, ifelse(rold_seed$variable  == "seed2012", 2012, ifelse(rold_seed$variable  == "seed2013", 2013, ifelse(rold_seed$variable  == "seed2014", 2014, ifelse(rold_seed$variable  == "seed2015", 2015, ifelse(rold_seed$variable  == "seed2016", 2016, NA)))))))))))

rold_seed1 <- rold_seed %>% 
  filter(!is.na(seed))
# View(rold_seed1)

rold_spike <- POAL_data_old_r %>% 
  rename("tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate("Birth Year" = case_when(Date == "Q2010" ~ "2010",
                                  Date != "Q2010" ~ Date))  %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike1_ = spike1, spike2_ = spike2) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", "spikelets_inflA5", "spikelets_infl13","spike1_","spike2_", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("spike1_",variable)~"A",
                              grepl("spike2_", variable)~"B"))
rold_spike$year<- ifelse(rold_spike$variable == "spike2007", 2007, ifelse(rold_spike$variable == "spike2008", 2008, ifelse(rold_spike$variable == "spike2009", 2009, ifelse(rold_spike$variable == "spike2010", 2010, ifelse(rold_spike$variable == "spike2011", 2011, ifelse(rold_spike$variable == "spikelets_inflA5", 2012, ifelse(rold_spike$variable  == "spikelets_infl13", 2013, ifelse(rold_spike$variable  == "spike1_", 2014, ifelse(rold_spike$variable == "spike2_", 2014, ifelse(rold_spike$variable  == "spikelets_infl15", 2015, ifelse(rold_spike$variable  == "spikelets_infl16", 2016, NA)))))))))))
rold_spike1 <- rold_spike %>% 
  filter(!is.na(spikelets))
# View(rold_spike1)

roldflw <- POAL_data_old_r %>%
  rename("tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate("Birth Year" = case_when(Date == "Q2010" ~ "2010",
                                  Date != "Q2010" ~ Date)) %>%   
  mutate(flw2007 = NA, flw2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "FLWTiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
roldflw$year<- ifelse(roldflw$variable == "flw2007", 2007, ifelse(roldflw$variable == "flw2008", 2008, ifelse(roldflw$variable == "FLWTiller09", 2009, ifelse(roldflw$variable == "FLWtiller10", 2010, ifelse(roldflw$variable == "FLWtiller11", 2011, ifelse(roldflw$variable  == "FLWtiller12", 2012, ifelse(roldflw$variable  == "FLWtiller13", 2013, ifelse(roldflw$variable  == "FLWtiller14", 2014, ifelse(roldflw$variable  == "FLWtiller15", 2015, ifelse(roldflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(roldflw)

rold_seedmerge_ss <- full_join(rold_spike1, rold_seed1, by = c( "plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year", "tillerid"))
# View(rold_seedmerge_ss)

rold_seedmerge_ssf <- merge(rold_seedmerge_ss, roldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
# View(rold_seedmerge_ssf)
# POAL(Old) recruits data includes tag # 10_19, which doesn't have any data collected for the reproductive data and is not present in the endo_demog_long file. This will be filtered out later when we merge the reproductive output with endo_demog_long 

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
# POALrepro <- POALrepro[!(is.na(POALrepro$flw)),]





# Combining seed productions measurements across years for the “New” POSY data -------------

## recoding for the year of measurement
## merging these measurements into one dataframe
po_seed <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA2", "seeds_InflB2", 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
         
po_seed$year<- ifelse(po_seed$variable == "seed2007", 2007, ifelse(po_seed$variable == "seed2008", 2008,ifelse(po_seed$variable == "seeds_InflA2", 2009, ifelse(po_seed$variable  == "seeds_InflB2", 2009, ifelse(po_seed$variable  == "seed2010", 2010, ifelse(po_seed$variable  == "seed2011", 2011, ifelse(po_seed$variable  == "seed2012", 2012, ifelse(po_seed$variable  == "seed2013", 2013,ifelse(po_seed$variable == "seed2014", 2014,ifelse(po_seed$variable == "seed2015", 2015,ifelse(po_seed$variable  == "seed2016", 2016, NA)))))))))))
po_seed1 <- po_seed %>% 
  filter(!is.na(seed))
# View(po_seed1)

po_spike <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spikelets_InflA2", "spikelets_InflB2", "spikelets_inflA3", 
                       "spikelets_inflB3", "spikelets_inflA4","spikelets_inflB4",
                       "spikelets_inflA5", "spikelets_inflA6", "spikelets_InflB6", 
                       "spikelets_inflC6", "spikelets_inflA7", "spikelets_InflB7", 
                       "spikelets_inflC7","spikelets_inflA8", "spikelets_InflB8", 
                       "spikelets_inflC8", "spike2016"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))

po_spike$year<- ifelse(po_spike$variable == "spike2007", 2007, ifelse(po_spike$variable == "spike2008", 2008,ifelse(po_spike$variable == "spikelets_InflA2", 2009, ifelse(po_spike$variable  == "spikelets_InflB2", 2009, ifelse(po_spike$variable  == "spikelets_inflA3", 2010, ifelse(po_spike$variable  == "spikelets_inflB3", 2010, ifelse(po_spike$variable  == "spikelets_inflB3", 2010, ifelse(po_spike$variable  == "spikelets_inflA4", 2011, ifelse(po_spike$variable == "spikelets_inflB4", 2011, ifelse(po_spike$variable == "spikelets_inflA5", 2012,ifelse(po_spike$variable == "spikelets_inflA6", 2013, ifelse(po_spike$variable == "spikelets_InflB6", 2013, ifelse(po_spike$variable == "spikelets_inflC6", 2013, ifelse(po_spike$variable  == "spikelets_inflA7", 2014, ifelse(po_spike$variable == "spikelets_InflB7", 2014, ifelse(po_spike$variable == "spikelets_inflC7", 2014, ifelse(po_spike$variable == "spikelets_inflA8", 2015, ifelse(po_spike$variable == "spikelets_InflB8", 2015, ifelse(po_spike$variable ==   "spikelets_inflC8", 2015, ifelse(po_spike$variable == "spike2016", 2016, NA))))))))))))))))))))
po_spike1 <- po_spike %>% 
  filter(!is.na(spikelets))
# View(po_spike1)


po_flw <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(flw2007 = 0, flw2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "flw2016"), 
       value.name = "flw") %>% 
  filter(!is.na(plot))
po_flw$year<- ifelse(po_flw$variable == "flw2007", 2007, ifelse(po_flw$variable == "Flwtillers1", 2008, ifelse(po_flw$variable  == "FlwTillers2", 2009, ifelse(po_flw$variable  == "FlwTillers3", 2010, ifelse(po_flw$variable  == "FlwTillers4", 2011, ifelse(po_flw$variable  == "FlwTillers5", 2012, ifelse(po_flw$variable  == "FlwTillers6", 2013,ifelse(po_flw$variable == "FlwTillers7", 2014,ifelse(po_flw$variable == "FlwTillers8", 2015,ifelse(po_flw$variable  == "flw2016", 2016, NA))))))))))
# View(po_flw)

po_seedmerge_ss <- full_join(po_spike1, po_seed1, by = c("plot", "pos", "tag", "Endo", 
                                                   "Birth Year","year","tillerid"))
# View(po_seedmerge_ss)


po_seedmerge_ssf <- merge(po_seedmerge_ss,po_flw, by = c("plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year"), all = TRUE)

# View(po_seedmerge_ssf)




# Combining repro measurements across years for the "New" POSY recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
po_rseed <- POSY_data_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2009", "seed2010", "seed2011",
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
po_rseed$year<- ifelse(po_rseed$variable == "seed2007", 2007, ifelse(po_rseed$variable == "seed2008", 2008, ifelse(po_rseed$variable == "seed2009", 2009, ifelse(po_rseed$variable  == "seed2010", 2010, ifelse(po_rseed$variable  == "seed2011", 2011, ifelse(po_rseed$variable  == "seed2012", 2012, ifelse(po_rseed$variable  == "seed2013", 2013,ifelse(po_rseed$variable == "seed2014", 2014,ifelse(po_rseed$variable == "seed2015", 2015,ifelse(po_rseed$variable  == "seed2016", 2016, NA))))))))))

po_rseed1 <- po_rseed %>% 
  filter(!is.na(seed))
# View(po_rseed1)

po_rspike <- POSY_data_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA,) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", 
                       "spikelets_inflA5", "spikelets_inflA5__1", 
                       "spikelets_infl14","spike1", "spike2", 
                       "spike3", "spike1__1", "spike2__1", "spike3__1"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                              grepl("spike2", variable) ~ "B",
                              grepl("spike3", variable) ~ "C"))
po_rspike$year<- ifelse(po_rspike$variable  == "spike2007", 2007, ifelse(po_rspike$variable  == "spike2008", 2008, ifelse(po_rspike$variable  == "spike2009", 2009, ifelse(po_rspike$variable  == "spike2010", 2010, ifelse(po_rspike$variable == "spike2011", 2011, ifelse(po_rspike$variable  == "spikelets_inflA5", 2012, ifelse(po_rspike$variable == "spikelets_inflA5__1", 2013, ifelse(po_rspike$variable == "spikelets_InflB6", 2013, ifelse(po_rspike$variable == "spikelets_inflC6", 2013, ifelse(po_rspike$variable  == "spikelets_infl14", 2014, ifelse(po_rspike$variable == "spike1", 2015, ifelse(po_rspike$variable == "spike2", 2015, ifelse(po_rspike$variable ==   "spike3", 2015, ifelse(po_rspike$variable == "spike1__1", 2016, ifelse(po_rspike$variable == "spike2__1", 2016, ifelse(po_rspike$variable == "spike3__1", 2016, NA))))))))))))))))

po_rspike1 <- po_rspike %>% 
  filter(!is.na(spikelets))# View(po_rspike)


po_rflw <- POSY_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = 0, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_rflw$year<- ifelse(po_rflw$variable == "flw2007", 2007, ifelse(po_rflw$variable == "flw2008", 2008, ifelse(po_rflw$variable == "flw2009", 2009, ifelse(po_rflw$variable == "FLWtiller10", 2010, ifelse(po_rflw$variable == "FLWtiller11", 2011, ifelse(po_rflw$variable  == "FLWtiller12", 2012, ifelse(po_rflw$variable  == "FLWtiller13", 2013, ifelse(po_rflw$variable  == "FLWtiller14", 2014, ifelse(po_rflw$variable  == "FLWtiller15", 2015, ifelse(po_rflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(po_rflw)


po_rseedmerge_ss <- full_join(po_rspike1, po_rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                           "Birth Year","year","tillerid"))
# View(po_rmerge_ss)

po_rseedmerge_ssf <- merge(po_rseedmerge_ss, po_rflw, by = c("plot", "pos", "tag", "Endo", 
                                                         "Birth Year","year"), all = TRUE)
# View(po_rseedmerge_ssf)






# Combining repro measurements across years for the “Old” POSY data ------------------


## Combining data across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe

po_oldseed <- POSY_data_old %>%
  mutate(seed2007 = NA, seed2008 = NA,  seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                 is.na(Date) ~ 2007)) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA3", "seeds_InflB3", "seed2010", 
                       "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

po_oldseed$year<- ifelse(po_oldseed$variable == "seed2007", 2007, ifelse(po_oldseed$variable == "seed2008", 2008, ifelse(po_oldseed$variable  == "seeds_InflA3", 2009, ifelse(po_oldseed$variable  == "seeds_InflB3",2009, ifelse(po_oldseed$variable == "seed2010", 2010, ifelse(po_oldseed$variable  == "seed2011", 2011, ifelse(po_oldseed$variable  == "seed2012", 2012, ifelse(po_oldseed$variable  == "seed2013", 2013,ifelse(po_oldseed$variable == "seed2014", 2014,ifelse(po_oldseed$variable == "seed2015", 2015,ifelse(po_oldseed$variable  == "seed2016", 2016, NA)))))))))))
po_oldseed1 <- po_oldseed %>% 
  filter(!is.na(seed))
# View(po_oldseed1)

po_oldspike <- POSY_data_old %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(spike2007 = NA, ) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("spike2007", "Spikelets_InflA1", "Spikelets_InflB1", 
                       "spikelets_InflA3","spikelets_InflB3", 
                       "spikelets_InflA4","spikelets_inflB4", "spikelets_InflA5", 
                       "spikelets_inflB5", "spikelets_inflA6", "spikelets_inflA7",
                       "spikelets_InflB7","spikelets_inflC7","spikelets_inflA8",
                       "spikelets_InflB8","spikelets_inflC8","spikelets_inflA9",
                       "spikelets_InflB9","spikelets_inflC9","spikelets_inflA10",
                       "spikelets_InflB10","spikelets_inflC10"), 
       value.name = "spikelets") %>% 
    mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                                grepl("B",variable) ~ "B",
                                grepl("C", variable) ~ "C"))

po_oldspike$year<- ifelse(po_oldspike$variable == "spike2007", 2007, ifelse(po_oldspike$variable == "Spikelets_InflA1", 2008,  ifelse(po_oldspike$variable == "Spikelets_InflB1", 2008, ifelse(po_oldspike$variable  == "spikelets_InflA3", 2009, ifelse(po_oldspike$variable  == "spikelets_InflB3", 2009, ifelse(po_oldspike$variable  == "spikelets_InflA4", 2010, ifelse(po_oldspike$variable  == "spikelets_inflB4", 2010, ifelse(po_oldspike$variable  == "spikelets_InflA5", 2011, ifelse(po_oldspike$variable  == "spikelets_inflB5", 2011, ifelse(po_oldspike$variable  == "spikelets_inflA6", 2012, ifelse(po_oldspike$variable  == "spikelets_inflA7", 2013, ifelse(po_oldspike$variable  == "spikelets_InflB7", 2013, ifelse(po_oldspike$variable  == "spikelets_inflC7", 2013, ifelse(po_oldspike$variable == "spikelets_inflA8", 2014, ifelse(po_oldspike$variable == "spikelets_InflB8", 2014, ifelse(po_oldspike$variable == "spikelets_inflC8", 2014, ifelse(po_oldspike$variable == "spikelets_inflA9", 2015,ifelse(po_oldspike$variable == "spikelets_InflB9", 2015, ifelse(po_oldspike$variable == "spikelets_inflC9", 2015, ifelse(po_oldspike$variable  == "spikelets_inflA10", 2016, ifelse(po_oldspike$variable  == "spikelets_InflB10", 2016, ifelse(po_oldspike$variable  == "spikelets_inflC10", 2016,NA))))))))))))))))))))))
po_oldspike1 <- po_oldspike %>% 
  filter(!is.na(spikelets), spikelets != 0)
# View(po_oldspike1)

po_oldflw <- POSY_data_old %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n",
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
po_oldflw$year<- ifelse(po_oldflw$variable == "flw2007", 2007, ifelse(po_oldflw$variable == "FlwTillers1", 2008, ifelse(po_oldflw$variable  == "FlwTillers3", 2009, ifelse(po_oldflw$variable  == "FlwTillers4", 2010, ifelse(po_oldflw$variable  == "FlwTillers5", 2011, ifelse(po_oldflw$variable  == "FlwTillers6", 2012, ifelse(po_oldflw$variable  == "FlwTillers7", 2013,ifelse(po_oldflw$variable == "FlwTillers8", 2014,ifelse(po_oldflw$variable == "FlwTillers9", 2015,ifelse(po_oldflw$variable  == "FlwTillers10", 2016, NA))))))))))
# View(po_oldflw)


po_oldseedmerge_ss <- full_join(po_oldseed1, po_oldspike1, by = c("plot", "pos", "tag", "Endo", 
                                                                "Birth Year","year","tillerid"))
# View(po_oldmerge_ss)

po_oldseedmerge_ssf <- merge(po_oldseedmerge_ss, po_oldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
# View(po_oldseedmerge_ssf)



# Combining  measurements across years for the “Old” POSY recruits data --------

## Combining data for recruits across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
po_roldseed <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("seed1", variable) ~ "A", 
                              grepl("seed2",variable) ~ "B",
                              grepl("seed3", variable) ~ "C"))
po_roldseed$year<- ifelse(po_roldseed$variable == "seed2007", 2007, ifelse(po_roldseed$variable == "seed2008", 2008, ifelse(po_roldseed$variable == "seed2009", 2009, ifelse(po_roldseed$variable == "seed2010", 2010, ifelse(po_roldseed$variable == "seed2011", 2011, ifelse(po_roldseed$variable  == "seed2012", 2012, ifelse(po_roldseed$variable  == "seed2013", 2013, ifelse(po_roldseed$variable  == "seed2014", 2014, ifelse(po_roldseed$variable  == "seed2015", 2015, ifelse(po_roldseed$variable  == "seed2016", 2016, NA))))))))))
po_roldseed1 <- po_roldseed %>% 
  filter(!is.na(seed))
# View(po_roldseed1)


po_roldspike <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 =NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", 
                       "spike2011", "spikelets_inflA12","spikelets_inflA13",
                       "spikelets_inflb13","spikelets_inflc13",
                       "spike1", "spike2", "spike3",
                       "spike1__1", "spike2__1", "spike3__1",
                       "spike1__2", "spike2__2", "spike3__2"),
       value.name = "spikelets") %>% 
    mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                                grepl("spike2",variable) ~ "B",
                                grepl("spike3", variable) ~ "C",
                                grepl("inflA", variable) ~ "A",
                                grepl("inflb", variable) ~ "B",
                                grepl("inflc", variable) ~ "C"))
  
po_roldspike$year<- ifelse(po_roldspike$variable == "spike2007", 2007, ifelse(po_roldspike$variable == "spike2008", 2008, ifelse(po_roldspike$variable == "spike2009", 2009, ifelse(po_roldspike$variable == "spike2010", 2010, ifelse(po_roldspike$variable == "spike2011", 2011, ifelse(po_roldspike$variable  == "spikelets_inflA12", 2012, ifelse(po_roldspike$variable  == "spikelets_inflA13", 2013, ifelse(po_roldspike$variable  == "spikelets_inflb13", 2013, ifelse(po_roldspike$variable  == "spikelets_inflc13", 2013, ifelse(po_roldspike$variable  == "spike1", 2014, ifelse(po_roldspike$variable  == "spike2", 2014, ifelse(po_roldspike$variable  == "spike3", 2014, ifelse(po_roldspike$variable  == "spike1__1", 2015, ifelse(po_roldspike$variable  == "spike2__1", 2015, ifelse(po_roldspike$variable  == "spike3__1", 2015, ifelse(po_roldspike$variable  == "spike1__2", 2016, ifelse(po_roldspike$variable  == "spike2__2", 2016,ifelse(po_roldspike$variable  == "spike3__2", 2016,NA))))))))))))))))))
po_roldspike1 <- po_roldspike %>% 
  filter(!is.na(spikelets))
#View(po_roldspike1)

po_roldflw <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWTiller12", 
                       "FLWTiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") %>% 
  distinct() # there are a couple of tag ids that are duplicated
po_roldflw$year<- ifelse(po_roldflw$variable == "flw2007", 2007, ifelse(po_roldflw$variable == "flw2008", 2008, ifelse(po_roldflw$variable == "flw2009", 2009, ifelse(po_roldflw$variable == "FLWtiller10", 2010, ifelse(po_roldflw$variable == "FLWtiller11", 2011, ifelse(po_roldflw$variable  == "FLWTiller12", 2012, ifelse(po_roldflw$variable  == "FLWTiller13", 2013, ifelse(po_roldflw$variable  == "FLWtiller14", 2014, ifelse(po_roldflw$variable  == "FLWtiller15", 2015, ifelse(po_roldflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(po_roldflw)

po_roldseedmerge_ss <-  full_join(po_roldspike1, po_roldseed1, by = c("plot", "pos", "tag", "Endo", 
                                                                   "Birth Year","year","tillerid"))
# View(po_roldseedmerge_ss)

po_roldseedmerge_ssf <- merge(po_roldseedmerge_ss, po_roldflw, by = c("plot", "pos", "tag", "Endo", 
                                                                     "Birth Year","year"), all = TRUE)
# View(po_roldseedmerge_ssf)






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

# POSYrepro <- POSYrepro[!(is.na(POSYrepro$flw)),]
# View(POSYrepro)





# Combining repro measurements across years for the LOAR data ------------------

## Combining measurements across years for reproduction measurements

# there is a separate sheet with the main raw data for seeds from 2008. It is pretty messy, but I am removing the columns that contain duplicated spikelet information as well as a few rows that have questionable data
LOAR_data_seed2008_long <- LOAR_data_seed2008 %>% 
  select(-contains("__1"), -contains("Avg"), -PropClaviceps1, -X__2, -X__3, -contains("Bug"), -tag2, -Notes,-SeedNotes1, -totunfilled) %>% 
  filter(!is.na(Tag), Spikelets != "-") %>% 
  group_by(Tag) %>% mutate(tillerid2008 = as.character(LETTERS[row_number()]),  year = 2008) %>% 
  mutate(seed2008 = as.numeric(`Unfilled Green`) + as.numeric(`Unfilled Brown`) + Filled) %>% 
  rename(spikelets2008 = Spikelets)

LOAR_data_seed2008_seed <- LOAR_data_seed2008_long %>% 
  select(Tag, tillerid2008, year, seed2008)
LOAR_data_seed2008_spike <- LOAR_data_seed2008_long %>% 
  select(Tag, tillerid2008, year, spikelets2008)
# there is a separate sheet with the main raw data for seeds from 2009. It is pretty messy, but I am removing the columns that contain duplicated spikelet information as well as a few rows that have questionable data
LOAR_data_seed2009_long <-LOAR_data_seed2009 %>% 
  select(-contains("__1"), -contains("Avg"), -contains("__2"), -contains("Bug"), -Notes, -TOTunfilled, -TOTUnfilled) %>% 
  filter(!is.na(Tag)) %>% 
  group_by(Tag) %>% mutate(tillerid2009 = as.character(LETTERS[row_number()]), year = 2009) %>% 
  mutate(seed2009 = `Unfilled Green` + `Unfilled Brown` + Filled) %>% 
  rename(spikelets2009 = Spikelets)

LOAR_data_seed2009_seed <- LOAR_data_seed2009_long %>% 
  select(Tag, tillerid2009, year, seed2009)
LOAR_data_seed2009_spike <- LOAR_data_seed2009_long %>% 
  select(Tag, tillerid2009, year, spikelets2009)

# I'm going to merge these two years into the spikelet and seed data respectively after producing the "long" version of the rest of the years.
# There is really no raw seed data in the LOAR excel sheet, just calculations from the spikelet counts, and some averages (some of which I think are taking the wrong cells, or omiting the *.12 rate of seeds per spikelet)
lseed <- LOAR_data%>% 
  rename("plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG" ) %>% 
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(seed2007 = NA,
         seed2008 = NA,
         seed2009 = NA,
         seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009",
                       "seed2010", "seed2011", "seed2012", "seed2013","seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
    mutate(tillerid = case_when(grepl("a", variable) ~ "A", 
                                grepl("b",variable) ~ "B",
                                grepl("c", variable) ~ "C"))
lseed$year <- ifelse(lseed$variable == "seed2007", 2007, ifelse(lseed$variable == "seed2008", 2008, ifelse(lseed$variable  == "seed2009", 2009, ifelse(lseed$variable  == "seed2010", 2010, ifelse(lseed$variable  == "seed2011", 2011, ifelse(lseed$variable  == "seed2012", 2012, ifelse(lseed$variable  == "seed2013", 2013,ifelse(lseed$variable == "seed2014", 2014,ifelse(lseed$variable == "seed2015", 2015,ifelse(lseed$variable  == "seed2016", 2016, NA))))))))))
# View(lseed)

# merging in the 2008 and 2009 seed data
lseed_1 <- full_join(lseed, LOAR_data_seed2008_seed, by = c("tag" = "Tag", "year" = "year"))


lseed_2 <- full_join(lseed_1, LOAR_data_seed2009_seed, by = c("tag" = "Tag", "year" = "year"))

lseed_2$seed_new <- ifelse(!is.na(lseed_2$seed), lseed_2$seed, ifelse(!is.na(lseed_2$seed2008), lseed_2$seed2008, ifelse(!is.na(lseed_2$seed2009), lseed_2$seed2009, NA)))
lseed_2$tillerid_new <- ifelse(!is.na(lseed_2$tillerid), lseed_2$tillerid, ifelse(!is.na(lseed_2$tillerid2008), lseed_2$tillerid2008, ifelse(!is.na(lseed_2$tillerid2009), lseed_2$tillerid2009, NA)))
lseed_3 <- lseed_2 %>% 
  select(plot, pos, tag, Endo, `Birth Year`, variable, year, seed_new, tillerid_new) %>% 
  rename(seed = seed_new, tillerid = tillerid_new) %>% 
  filter(!is.na(seed))



lspike <- LOAR_data %>% 
  rename("plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(spike2007 = NA,
         spike2008 = NA,
         spike2009 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("spike2007", "spike2008", "spike2009",
                       "spikelets_inflA3", "spikelets_inflB3",
                       "spikelets_inflA5", "spikelets_inflB5",
                       "spikelets_inflA6", "spikelets_inflB6",
                       "spikelets_inflA7", "spikelets_inflB7",
                       "spikelets_inflA8", "spikelets_inflB8",
                       "spikelets_inflA9", "spikelets_inflB9"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
lspike$year<- ifelse(lspike$variable == "spike2007", 2007, ifelse(lspike$variable == "spike2008", 2008,  ifelse(lspike$variable == "spike2009", 2009, ifelse(lspike$variable  == "spikelets_inflA3", 2010, ifelse(lspike$variable  == "spikelets_inflB3", 2010, ifelse(lspike$variable  == "spikelets_inflA5", 2012, ifelse(lspike$variable  == "spikelets_inflB5", 2012, ifelse(lspike$variable  == "spikelets_inflA6", 2013, ifelse(lspike$variable  == "spikelets_inflB6", 2013, ifelse(lspike$variable  == "spikelets_inflA7", 2014, ifelse(lspike$variable  == "spikelets_inflB7", 2014, ifelse(lspike$variable == "spikelets_inflA8", 2015, ifelse(lspike$variable == "spikelets_inflB8", 2015, ifelse(lspike$variable == "spikelets_inflA9", 2016, ifelse(lspike$variable == "spikelets_inflB9", 2016, NA)))))))))))))))

# View(lspike)

lspike_1 <- full_join(lspike, LOAR_data_seed2008_spike, by = c("tag" = "Tag", "year" = "year"))
lspike_2 <- full_join(lspike_1, LOAR_data_seed2009_spike, by = c("tag" = "Tag", "year" = "year"))

lspike_2$spikelets_new <- ifelse(!is.na(lspike_2$spikelets), lspike_2$spikelets, ifelse(!is.na(lspike_2$spikelets2008), lspike_2$spikelets2008, ifelse(!is.na(lspike_2$spikelets2009), lspike_2$spikelets2009, NA)))
lspike_2$tillerid_new <- ifelse(!is.na(lspike_2$tillerid), lspike_2$tillerid, ifelse(!is.na(lspike_2$tillerid2008), lspike_2$tillerid2008, ifelse(!is.na(lspike_2$tillerid2009), lspike_2$tillerid2009, NA)))
lspike_3 <- lspike_2 %>% 
  select(plot, pos, tag, Endo, `Birth Year`, variable, year, spikelets_new, tillerid_new) %>% 
  rename(spikelets = spikelets_new, tillerid = tillerid_new) %>% 
  filter(!is.na(spikelets))


lflw <- LOAR_data %>% 
  rename("plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "FlwTillers1", "FLWTiller2", "FLWTillers3", "FLWTillers4",
                       "FLWTillers5", "FLWTillers6", "FLWTillers7", 
                       "FLWTillers8", "FLWTillers9"), 
       value.name = "flw") 
lflw$year<- ifelse(lflw$variable == "flw2007", 2007, ifelse(lflw$variable == "FlwTillers1", 2008, ifelse(lflw$variable  == "FLWTiller2", 2009, ifelse(lflw$variable  == "FLWTillers3", 2010, ifelse(lflw$variable  == "FLWTillers4", 2011, ifelse(lflw$variable  == "FLWTillers5", 2012, ifelse(lflw$variable  == "FLWTillers6", 2013,ifelse(lflw$variable == "FLWTillers7", 2014,ifelse(lflw$variable == "FLWTillers8", 2015,ifelse(lflw$variable  == "FLWTillers9", 2016, NA))))))))))
# View(lflw)

l_seedmerge_ss <- full_join(lspike_3, lseed_3, by = c("plot", "pos", "tag", "Endo", 
                                                      "Birth Year", "year", "tillerid"))
# View(l_seedmerge_ss)


l_seedmerge_ssf <- merge(l_seedmerge_ss, lflw, by = c("plot", "pos", "tag", "Endo", 
                                                      "Birth Year", "year"), all = TRUE)

# View(l_seedmerge_ssf)





# Combining repro measurements across years for the LOAR recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
l_rseed <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID", "Endo" = "endo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed")   %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
l_rseed$year<- ifelse(l_rseed$variable == "seed2007", 2007, ifelse(l_rseed$variable == "seed2008", 2008, ifelse(l_rseed$variable == "seed2009", 2009,ifelse(l_rseed$variable == "seed2010", 2010, ifelse(l_rseed$variable == "seed2011", 2011, ifelse(l_rseed$variable  == "seed2012", 2012, ifelse(l_rseed$variable  == "seed2013", 2013, ifelse(l_rseed$variable  == "seed2014", 2014, ifelse(l_rseed$variable  == "seed2015", 2015, ifelse(l_rseed$variable  == "seed2016", 2016, NA))))))))))
l_rseed_1 <- l_rseed %>% 
  filter(!is.na(seed))
# View(l_rseed_1)


l_rspike <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID", "Endo" = "endo") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", 
                       "spike2010","spike2011", "spike2012",
                       "spike1", "spike2", "spike1__1", "spike2__1",
                       "spike1_15", "spike2_15", "spike1_16", "spike2_16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                              grepl("spike2",variable) ~ "B"))
l_rspike$year<- ifelse(l_rspike$variable == "spike2007", 2007, ifelse(l_rspike$variable == "spike2008", 2008, ifelse(l_rspike$variable == "spike2009", 2009, ifelse(l_rspike$variable == "spike2010", 2010, ifelse(l_rspike$variable == "spike2011", 2011, ifelse(l_rspike$variable  == "spike2012", 2012, ifelse(l_rspike$variable  == "spike1", 2013, ifelse(l_rspike$variable  == "spike2", 2013, ifelse(l_rspike$variable  == "spike1__1", 2014, ifelse(l_rspike$variable  == "spike2__1", 2014, ifelse(l_rspike$variable  == "spike1_15", 2015,ifelse(l_rspike$variable  == "spike2_15", 2015, ifelse(l_rspike$variable  == "spike1_16", 2016,ifelse(l_rspike$variable  == "spike2_16", 2016, NA))))))))))))))
l_rspike_1 <- l_rspike %>% 
  filter(!is.na(spikelets))
# View(l_rspike_1)

l_rflw <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID", "Endo" = "endo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FlwTillers10","Flw11", "Flw12", 
                       "Flw13", "Flw14", "Flw15",
                       "Flw16"),
       value.name = "flw") 
l_rflw$year<- ifelse(l_rflw$variable == "flw2007", 2007, ifelse(l_rflw$variable == "flw2008", 2008, ifelse(l_rflw$variable == "flw2009", 2009, ifelse(l_rflw$variable == "FlwTillers10", 2010, ifelse(l_rflw$variable == "Flw11", 2011, ifelse(l_rflw$variable  == "Flw12", 2012, ifelse(l_rflw$variable  == "Flw13", 2013, ifelse(l_rflw$variable  == "Flw14", 2014, ifelse(l_rflw$variable  == "Flw15", 2015, ifelse(l_rflw$variable  == "Flw16", 2016, NA))))))))))
# View(l_rflw)

l_rseedmerge_ss <- full_join(l_rspike_1, l_rseed_1, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Birth Year", "year", "tillerid"))
# View(l_rseedmerge_ss)

l_rseedmerge_ssf <- merge(l_rseedmerge_ss, l_rflw, by = c( "plot", "pos", "tag", "Endo", 
                                                    "Birth Year", "year"),all = TRUE)
# View(l_rseedmerge_ssf)

# Combining recruit and original plant data
l_seedmerge_ssf <- l_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                       "tillerid", "flw","spikelets", "seed")]
l_rseedmerge_ssf <- l_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                         "tillerid", "flw","spikelets", "seed")]


LOARrepro <- l_seedmerge_ssf %>% 
  rbind(l_rseedmerge_ssf) %>% 
  mutate(species = "LOAR")

# LOARrepro <- LOARrepro[!(is.na(LOARrepro$flw)),]
# View(LOARrepro)










# Combining repro measurements across years for the FESU data ------------------

## Combining measurements across years for Seed, Spike, and Flowering using melt
## Recoding those measurements for the year they are taken

fseed <- FESU_data %>%
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA2", "seeds_InflB2",  
                       "seed2010", "seed2011","seed2012",
                       "seed2013", "seed2014", "seed2015", 
                       "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

fseed$year<- ifelse(fseed$variable == "seed2007", 2007, ifelse(fseed$variable == "seed2008", 2008, ifelse(fseed$variable  == "seeds_InflA2", 2009, ifelse(fseed$variable  == "seeds_InflB2", 2009, ifelse(fseed$variable  == "seed2010", 2010, ifelse(fseed$variable  == "seed2011", 2011, ifelse(fseed$variable  == "seed2012", 2012, ifelse(fseed$variable  == "seed2013", 2013,ifelse(fseed$variable == "seed2014", 2014,ifelse(fseed$variable == "seed2015", 2015,ifelse(fseed$variable  == "seed2016", 2016, NA)))))))))))
fseed1 <- fseed %>% 
  filter(!is.na(seed))
# View(fseed1)

fspike <- FESU_data %>% 
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(spike2007 = NA, spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "spike2008", "spikelets_inflA2", "spikelets_InflB2",
                       "spikelets_inflA3", "spikelets_inflB3",
                       "spikelets_inflA4", "spikelets_inflB4",
                       "spikelets_inflA5", "spikelets_inflB5",
                       "spikelets_inflA6", "spikelets_inflB6",
                       "spikelets_inflA7", "spikelets_inflB7",
                       "spikelets_inflA8", "spikelets_inflB8",
                       "spikelets_inflA9", "spikelets_inflB9"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

fspike$year<- ifelse(fspike$variable == "spike2007", 2007, ifelse(fspike$variable == "spike2008", 2008, ifelse(fspike$variable  == "spikelets_inflA2", 2009, ifelse(fspike$variable  == "spikelets_InflB2", 2009, ifelse(fspike$variable  == "spikelets_inflA3", 2010, ifelse(fspike$variable  == "spikelets_inflB3", 2010, ifelse(fspike$variable  == "spikelets_inflA4", 2011, ifelse(fspike$variable  == "spikelets_inflB4", 2011, ifelse(fspike$variable  == "spikelets_inflA5", 2012, ifelse(fspike$variable  == "spikelets_inflB5", 2012, ifelse(fspike$variable  == "spikelets_inflA6", 2013, ifelse(fspike$variable  == "spikelets_inflB6", 2013, ifelse(fspike$variable == "spikelets_inflA7", 2014, ifelse(fspike$variable == "spikelets_inflB7", 2014, ifelse(fspike$variable == "spikelets_inflA8", 2015, ifelse(fspike$variable == "spikelets_inflB8", 2015, ifelse(fspike$variable  == "spikelets_inflA9", 2016, ifelse(fspike$variable  == "spikelets_inflB9", 2016, NA))))))))))))))))))
fspike1 <- fspike %>% 
  filter(!is.na(spikelets))
# View(fspike1)


fflw <- FESU_data %>% 
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = 0, flw2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "flw2008", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
fflw$year<- ifelse(fflw$variable == "flw2007", 2007, ifelse(fflw$variable == "flw2008", 2008, ifelse(fflw$variable  == "FlwTillers2", 2009, ifelse(fflw$variable  == "FlwTillers3", 2010, ifelse(fflw$variable  == "FlwTillers4", 2011, ifelse(fflw$variable  == "FlwTillers5", 2012, ifelse(fflw$variable  == "FlwTillers6", 2013,ifelse(fflw$variable == "FlwTillers7", 2014,ifelse(fflw$variable == "FlwTillers8", 2015,ifelse(fflw$variable  == "FlwTillers9", 2016, NA))))))))))
# View(fflw)

f_seedmerge_ss <- full_join(fspike1, fseed1, by = c("plot", "pos", "tag", "Endo", 
                                                  "Birth Year","year","tillerid"))
# View(f_seedmerge_ss)

f_seedmerge_ssf <- merge(f_seedmerge_ss, fflw, by = c("plot", "pos", "tag", "Endo", 
                                                      "Birth Year","year"), all = TRUE)
# View(f_seedmerge_ssf)


# Combining repro measurements across years for the FESU recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
f_rseed <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

f_rseed$year<- ifelse(f_rseed$variable == "seed2007", 2007, ifelse(f_rseed$variable == "seed2008", 2008, ifelse(f_rseed$variable == "seed2009", 2009, ifelse(f_rseed$variable == "seed2010", 2010, ifelse(f_rseed$variable == "seed2011", 2011, ifelse(f_rseed$variable  == "seed2012", 2012, ifelse(f_rseed$variable  == "seed2013", 2013, ifelse(f_rseed$variable  == "seed2014", 2014, ifelse(f_rseed$variable  == "seed2015", 2015, ifelse(f_rseed$variable  == "seed2016", 2016, NA))))))))))
f_rseed1 <- rseed %>% 
  filter(!is.na(seed))
# View(f_rseed1)

f_rspike <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012",
                       "spikelets_infla13", "spikelets_inflb13",
                       "spikelets_infla14", "spikelets_inflb14",
                       "spikelets_infla15", "spikelets_inflb15",
                       "spikelets_infla16", "spikelets_inflb16",
                       "spikelets_inflc16"
                       ),
       value.name = "spikelets")   %>% 
  mutate(tillerid = case_when(grepl("a", variable) ~ "A", 
                              grepl("b",variable) ~ "B",
                              grepl("c", variable) ~ "C"))

f_rspike$year<- ifelse(f_rspike$variable == "spike2007", 2007, ifelse(f_rspike$variable == "spike2008", 2008, ifelse(f_rspike$variable == "spike2009", 2009, ifelse(f_rspike$variable == "spike2010", 2010, ifelse(f_rspike$variable == "spike2011", 2011, ifelse(f_rspike$variable  == "spike2012", 2012, ifelse(f_rspike$variable  == "spikelets_infla13", 2013, ifelse(f_rspike$variable  == "spikelets_inflb13", 2013, ifelse(f_rspike$variable  == "spikelets_infla14", 2014, ifelse(f_rspike$variable  == "spikelets_inflb14", 2014, ifelse(f_rspike$variable  == "spikelets_infla15", 2015, ifelse(f_rspike$variable  == "spikelets_inflb15", 2015, ifelse(f_rspike$variable  == "spikelets_infla16", 2016, ifelse(f_rspike$variable  == "spikelets_inflb16", 2016, ifelse(f_rspike$variable  == "spikelets_inflc16", 2016,NA)))))))))))))))
f_rspike1 <- f_rspike %>% 
  filter(!is.na(spikelets))
# View(f_rspike1)

f_rflw <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLW10", "FLW11", "FLW12", 
                       "FLW13", "FLW14", "FLW15",
                       "FLW16"),
       value.name = "flw") %>% 
  distinct() # There are a few duplicated rows of tag id's
f_rflw$year<- ifelse(f_rflw$variable == "flw2007", 2007, ifelse(f_rflw$variable == "flw2008", 2008,ifelse(f_rflw$variable == "flw2009", 2009,ifelse(f_rflw$variable == "FLW10", 2010, ifelse(f_rflw$variable == "FLW11", 2011, ifelse(f_rflw$variable  == "FLW12", 2012, ifelse(f_rflw$variable  == "FLW13", 2013, ifelse(f_rflw$variable  == "FLW14", 2014, ifelse(f_rflw$variable  == "FLW15", 2015, ifelse(f_rflw$variable  == "FLW16", 2016, NA))))))))))
# View(f_rflw)

f_rseedmerge_ss <- full_join(f_rspike1, f_rseed1, by = c("plot", "pos", "tag", "Endo", 
                                                       "Birth Year","year","tillerid"))
# View(f_rseedmerge_ss)

f_rseedmerge_ssf <- merge(f_rseedmerge_ss, f_rflw, by = c("plot", "pos", "tag", "Endo", 
                                                          "Birth Year","year"), all = TRUE)
# View(f_rseedmerge_ssf)








# Combining the original and recruit FESU repro dataframes ----------------------
f_seedmerge_ssf <- f_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                             "tillerid", "flw","spikelets", "seed")]

f_rseedmerge_ssf <- f_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                               "tillerid", "flw","spikelets", "seed")]

FESUrepro <- f_seedmerge_ssf %>% 
  rbind(f_rseedmerge_ssf) %>% 
  mutate(species = "FESU")

# FESUrepro <- FESUrepro[!(is.na(FESUrepro$flw)),]

# View(FESUrepro)






# Combining repro measurements across years for the ELVI data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken
ELVI_seed_tiller <- merge(ELVI_data, ELVI_data_seed, by = "TAG") #there is a separate sheet with the main raw data for seeds.

elviseed <- ELVI_seed_tiller %>%
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(seed2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "Seeds_Infl08", 
                       "INfl1Seeds09","Infl2Seeds09","Infl3Seeds09","Infl4Seeds09", 
                       "Infl1Seeds10","Infl2Seeds10",
                       "Infl1Seeds11","Infl2Seeds11",
                       "Seeds1_Infl12", "Seeds2_Infl12",
                       "Seeds1_Infl13", "Seeds2_Infl13",
                       "Seeds1_Infl14", "Seeds2_Infl14",
                       "Seeds1_Infl15", "Seeds2_Infl15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "seed")    %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1S", variable) ~ "A", 
                              grepl("fl2S",variable) ~ "B",
                              grepl("fl3S", variable) ~ "C",
                              grepl("fl4S", variable) ~ "D",
                              grepl("Seeds_", variable) ~ "A",
                              grepl("Seeds1", variable) ~ "A",
                              grepl("Seeds2", variable) ~ "B"))

elviseed$year<- ifelse(elviseed$variable == "seed2007", 2007, ifelse(elviseed$variable == "Seeds_Infl08", 2008, ifelse(elviseed$variable  == "INfl1Seeds09", 2009, ifelse(elviseed$variable  == "Infl2Seeds09", 2009, ifelse(elviseed$variable  == "Infl3Seeds09", 2009, ifelse(elviseed$variable  == "Infl4Seeds09", 2009, ifelse(elviseed$variable  == "Infl1Seeds10", 2010, ifelse(elviseed$variable  == "Infl2Seeds10", 2010,ifelse(elviseed$variable  == "Infl1Seeds11", 2011,ifelse(elviseed$variable  == "Infl2Seeds11", 2011, ifelse(elviseed$variable  ==  "Seeds1_Infl12", 2012, ifelse(elviseed$variable  ==  "Seeds2_Infl12", 2012, ifelse(elviseed$variable  ==  "Seeds1_Infl13", 2013, ifelse(elviseed$variable  ==  "Seeds2_Infl13", 2013, ifelse(elviseed$variable =="Seeds1_Infl14", 2014, ifelse(elviseed$variable =="Seeds2_Infl14", 2014, ifelse(elviseed$variable == "Seeds1_Infl15", 2015, ifelse(elviseed$variable == "Seeds2_Infl15", 2015, ifelse(elviseed$variable  =="Seeds1_Infl16", 2016, ifelse(elviseed$variable  =="Seeds2_Infl16", 2016, NA))))))))))))))))))))
elviseed1 <- elviseed %>% 
  filter(!is.na(seed))
# View(elviseed1)

elvispike <- ELVI_seed_tiller %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012","spike2013", 
                       "spike2014", "spike2015", "spike2016"), 
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))

elvispike$year<- ifelse(elvispike$variable == "spike2007", 2007, ifelse(elvispike$variable == "spike2008", 2008, ifelse(elvispike$variable  == "spike2009", 2009, ifelse(elvispike$variable  == "spike2010", 2010, ifelse(elvispike$variable  == "spike2011", 2011, ifelse(elvispike$variable  == "spike2012", 2012, ifelse(elvispike$variable  == "spike2013", 2013, ifelse(elvispike$variable == "spike2014", 2014, ifelse(elvispike$variable == "spike2015", 2015, ifelse(elvispike$variable  == "spike2016", 2016, NA))))))))))
elvispike1 <- elvispike %>% 
  filter(!is.na(spikelets))
# View(elvispike1)

elviflw <- ELVI_seed_tiller %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = 0) %>% 
 melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers08", "FlwTillers09", "FlwTillers10", 
                       "FlwTillers11", "FlwTillers12", "FlwTillers13", 
                       "FlwTillers14", "FlwTillers15", "FlwTillers16"), 
       value.name = "flw") 
elviflw$year<- ifelse(elviflw$variable == "flw2007", 2007,ifelse(elviflw$variable == "FlwTillers08", 2008, ifelse(elviflw$variable  == "FlwTillers09", 2009, ifelse(elviflw$variable  == "FlwTillers10", 2010, ifelse(elviflw$variable  == "FlwTillers11", 2011, ifelse(elviflw$variable  == "FlwTillers12", 2012, ifelse(elviflw$variable  == "FlwTillers13", 2013,ifelse(elviflw$variable == "FlwTillers14", 2014,ifelse(elviflw$variable == "FlwTillers15", 2015,ifelse(elviflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elviflw)

elvi_seedmerge_ss <- full_join(elvispike1, elviseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Birth Year","year","tillerid"))
# View(elvi_seedmerge_ss)

elvi_seedmerge_ssf <- merge(elvi_seedmerge_ss, elviflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
# View(elvi_seedmerge_ssf)



# Combining repro measurements across years for the ELVI recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
elvi_rseed <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011", "seed2012", 
                       "seed2013", "seed2014", 
                       "seeds/infl1/15","seeds/infl2/15", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B"))

elvi_rseed$year<- ifelse(elvi_rseed$variable == "seed2007", 2007,ifelse(elvi_rseed$variable == "seed2008", 2008,ifelse(elvi_rseed$variable == "seed2009", 2009, ifelse(elvi_rseed$variable == "seed2010", 2010, ifelse(elvi_rseed$variable == "seed2011", 2011, ifelse(elvi_rseed$variable  == "seed2012", 2012, ifelse(elvi_rseed$variable  == "seed2013", 2013, ifelse(elvi_rseed$variable  == "seed2014", 2014, ifelse(elvi_rseed$variable  == "seeds/infl1/15", 2015, ifelse(elvi_rseed$variable  == "seeds/infl2/15", 2015, ifelse(elvi_rseed$variable == "seed2016", 2016, NA)))))))))))
elvi_rseed1 <- elvi_rseed %>% 
  filter(!is.na(seed))
# View(elvi_rseed1)

elvi_rspike <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010","spike2011", "spike2012", 
                       "spike2013", "spike2014", "spike2015", "spike2016"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))
elvi_rspike$year<- ifelse(elvi_rspike$variable == "spike2007", 2007, ifelse(elvi_rspike$variable == "spike2008", 2008, ifelse(elvi_rspike$variable == "spike2009", 2009, ifelse(elvi_rspike$variable == "spike2010", 2010, ifelse(elvi_rspike$variable == "spike2011", 2011, ifelse(elvi_rspike$variable  == "spike2012", 2012, ifelse(elvi_rspike$variable  == "spike2013", 2013, ifelse(elvi_rspike$variable  == "spike2014", 2014, ifelse(elvi_rspike$variable  == "spike2015", 2015, NA)))))))))
elvi_rspike1 <- elvi_rspike %>% 
  filter(!is.na(spikelets))
# View(elvi_rspike1)

elvi_rflw <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = NA, flw2008 =NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "FlwTillers09", "FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15"),
       value.name = "flw") 
elvi_rflw$year<- ifelse(elvi_rflw$variable == "flw2007", 2007, ifelse(elvi_rflw$variable == "flw2008", 2008, ifelse(elvi_rflw$variable == "FlwTillers09", 2009, ifelse(elvi_rflw$variable == "FlwTillers10", 2010, ifelse(elvi_rflw$variable == "FlwTillers11", 2011, ifelse(elvi_rflw$variable  == "FlwTillers12", 2012, ifelse(elvi_rflw$variable  == "FlwTillers13", 2013, ifelse(elvi_rflw$variable  == "FlwTillers14", 2014, ifelse(elvi_rflw$variable  == "FlwTillers15", 2015, ifelse(elvi_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elvi_rflw)

elvi_rseedmerge_ss <- full_join(elvi_rspike1, elvi_rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                                 "Birth Year","year","tillerid"))
# View(elvi_rseedmerge_ss)

elvi_rseedmerge_ssf <- merge(elvi_rseedmerge_ss, elvi_rflw, by = c("plot", "pos", "tag", "Endo", 
                                                                   "Birth Year","year"), all = TRUE)
# View(elvi_rseedmerge_ssf)




# Combining the  original and recruit ELVI repro dataframes ---------
elvi_seedmerge_ssf <- elvi_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                           "tillerid", "flw","spikelets", "seed")]

elvi_rseedmerge_ssf <- elvi_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                             "tillerid", "flw","spikelets", "seed")]


ELVIrepro <- elvi_seedmerge_ssf %>% 
  rbind(elvi_rseedmerge_ssf) %>% 
  mutate(species = "ELVI")

# ELVIrepro <- ELVIrepro[!(is.na(ELVIrepro$flw)),]

# View(ELVIrepro)








# Combining repro measurements across years for the ELRI data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken
ELRI_seed_tiller <- merge(ELRI_data, ELRI_data_seed2009, by.x = c("PLOT", "POS", "TAG", "ENDO", "Planted Date","TRT", "Loc'n", "Plant"), by.y = c("plot", "pos", "tag", "Endo", "Planted Date","TRT", "Loc'n", "Plant"), all.x = TRUE) #there is a separate sheet with the main raw data for seeds.


elriseed <- ELRI_seed_tiller %>%
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>% 
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(seed2007 = NA, seed2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", 
                       "Infl1Seeds2", "Infl2Seeds2", "Infl3Seeds2", "Infl4Seeds2", 
                       "Seeds_Infl10","seeds_1_11","seeds_2_11",
                       "seeds_1_12","seeds_2_12",
                       "seeds_1_13","seeds_2_13",
                       "seeds_1_14","seeds_2_14",
                       "seeds_1_15","seeds_2_15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("_1_", variable) ~ "A",
                              grepl("_2_", variable) ~ "B", 
                              grepl("fl1S",variable) ~ "A",
                              grepl("fl2S",variable)~ "B",
                              grepl("fl3S", variable) ~ "C",
                              grepl("fl4S", variable) ~ "D",
                              grepl("Infl10", variable) ~ "A"))

elriseed$year<- ifelse(elriseed$variable == "seed2007", 2007, ifelse(elriseed$variable == "seed2008", 2008, ifelse(elriseed$variable  == "Infl1Seeds2", 2009, ifelse(elriseed$variable  == "Infl2Seeds2", 2009, ifelse(elriseed$variable  == "Infl3Seeds2", 2009, ifelse(elriseed$variable  == "Infl4Seeds2", 2009, ifelse(elriseed$variable  == "Seeds_Infl10", 2010, ifelse(elriseed$variable  == "seeds_1_11", 2011, ifelse(elriseed$variable  == "seeds_2_11", 2011, ifelse(elriseed$variable  == "seeds_1_12", 2012, ifelse(elriseed$variable  == "seeds_2_12", 2012,ifelse(elriseed$variable  == "seeds_1_13", 2013,ifelse(elriseed$variable  == "seeds_2_13", 2013, ifelse(elriseed$variable == "seeds_1_14", 2014,ifelse(elriseed$variable == "seeds_2_14", 2014, ifelse(elriseed$variable == "seeds_1_15", 2015, ifelse(elriseed$variable == "seeds_2_15", 2015,ifelse(elriseed$variable  == "Seeds1_Infl16", 2016,ifelse(elriseed$variable  == "Seeds2_Infl16", 2016, NA)))))))))))))))))))
elriseed1 <- elriseed %>% 
  filter(!is.na(seed), seed != ".")
# View(elriseed1)

elrispike <- ELRI_seed_tiller %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT"), 
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012","spike2013", 
                       "spike2014", "spike2015", "spike2016"), 
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("_1_", variable) ~ "A",
                              grepl("_2_", variable) ~ "B", 
                              grepl("fl1",variable) ~ "A",
                              grepl("fl2",variable)~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))
elrispike$year<- ifelse(elrispike$variable == "spike2007", 2007, ifelse(elrispike$variable == "spike2008", 2008, ifelse(elrispike$variable  == "spike2009", 2009, ifelse(elrispike$variable  == "spike2010", 2010, ifelse(elrispike$variable  == "spike2011", 2011, ifelse(elrispike$variable  == "spike2012", 2012, ifelse(elrispike$variable  == "spike2013", 2013, ifelse(elrispike$variable == "spike2014", 2014, ifelse(elrispike$variable == "spike2015", 2015, ifelse(elrispike$variable  == "spike2016", 2016, NA))))))))))
elrispike1 <- elrispike %>% 
  filter(!is.na(spikelets))
# View(elrispike1)

elriflw <- ELRI_seed_tiller %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers08", "FlwTillers09.x", "FLwTillers10", 
                       "FlwTiller11", "FlwTiller12", "FlwTiller13", 
                       "FlwTiller14", "FlwTiller15", "FlwTillers16"), 
       value.name = "flw") 
elriflw$year<- ifelse(elriflw$variable == "flw2007", 2007, ifelse(elriflw$variable == "FlwTillers08", 2008, ifelse(elriflw$variable  == "FlwTillers09.x", 2009, ifelse(elriflw$variable  == "FLwTillers10", 2010, ifelse(elriflw$variable  == "FlwTiller11", 2011, ifelse(elriflw$variable  == "FlwTiller12", 2012, ifelse(elriflw$variable  == "FlwTiller13", 2013,ifelse(elriflw$variable == "FlwTiller14", 2014,ifelse(elriflw$variable == "FlwTiller15", 2015,ifelse(elriflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elriflw)

elri_seedmerge_ss <- full_join(elrispike1, elriseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                            "Birth Year","year","tillerid"))
# View(elri_seedmerge_ss)

elri_seedmerge_ssf <- merge(elri_seedmerge_ss, elriflw, by = c( "plot", "pos", "tag", "Endo", 
                                                                "Birth Year","year"), all = TRUE)
# View(elri_merge_ssf)


# Combining repro measurements across years for the ELRI recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
elri_rseed <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", 
                       "seed2011", "seed2012", "seed2013",
                       "seeds1_14","seeds2_14", "seeds3_14",
                       "seeds1_15","seeds2_15", "seeds3_15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("seeds1_", variable) ~ "A",
                              grepl("seeds2_", variable) ~ "B"))
elri_rseed$year<- ifelse(elri_rseed$variable == "seed2007", 2007,ifelse(elri_rseed$variable == "seed2008", 2008,ifelse(elri_rseed$variable == "seed2009", 2009,ifelse(elri_rseed$variable == "seed2010", 2010, ifelse(elri_rseed$variable == "seed2011", 2011, ifelse(elri_rseed$variable  == "seed2012", 2012, ifelse(elri_rseed$variable  == "seed2013", 2013, ifelse(elri_rseed$variable  == "seeds1_14", 2014, ifelse(elri_rseed$variable  == "seeds2_14", 2014, ifelse(elri_rseed$variable  == "seeds3_14", 2014, ifelse(elri_rseed$variable  == "seeds1_15", 2015, ifelse(elri_rseed$variable  == "seeds2_15", 2015, ifelse(elri_rseed$variable  == "seeds3_15", 2015, ifelse(elri_rseed$variable == "Seeds1_Infl16", 2016, ifelse(elri_rseed$variable == "Seeds2_Infl16", 2016,NA)))))))))))))))
elri_rseed1 <- elri_rseed %>% 
  filter(!is.na(seed))
# View(elri_rseed1)

elri_rspike <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012","spikelets_infla13","spikelets_inflb13", 
                       "spike2014", "spike2015", "spike2016"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("seeds1_", variable) ~ "A",
                              grepl("seeds2_", variable) ~ "B",
                              grepl("infla", variable) ~ "A",
                              grepl("inflb", variable) ~ "B"))
elri_rspike$year<- ifelse(elri_rspike$variable == "seed2007", 2007, ifelse(elri_rspike$variable == "seed2008", 2008, ifelse(elri_rspike$variable == "seed2009", 2009, ifelse(elri_rspike$variable == "seed2010", 2010, ifelse(elri_rspike$variable == "seed2011", 2011, ifelse(elri_rspike$variable  == "seed2012", 2012, ifelse(elri_rspike$variable  == "spikelets_infla13", 2013, ifelse(elri_rspike$variable  == "spikelets_inflb13", 2013, ifelse(elri_rspike$variable  == "seed2014", 2014, ifelse(elri_rspike$variable  == "seed2015", 2015, ifelse(elri_rspike$variable == "seed2016", 2016, NA)))))))))))
elri_rspike1 <- elri_rspike %>% 
  filter(!is.na(spikelets))
# View(elri_rspike1)

elri_rflw <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "FLWtiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLW13", "FLW14", "FLW15", "FlwTillers16"),
       value.name = "flw") 
elri_rflw$year<- ifelse(elri_rflw$variable == "seed2007", 2007, ifelse(elri_rflw$variable == "seed2008", 2008, ifelse(elri_rflw$variable == "FLWtiller09", 2009, ifelse(elri_rflw$variable == "FLWtiller10", 2010, ifelse(elri_rflw$variable == "FLWtiller11", 2011, ifelse(elri_rflw$variable  == "FLWtiller12", 2012, ifelse(elri_rflw$variable  == "FLW13", 2013, ifelse(elri_rflw$variable  == "FLW14", 2014, ifelse(elri_rflw$variable  == "FLW15", 2015, ifelse(elri_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elri_rflw)

elri_rseedmerge_ss <- full_join(elri_rspike1, elri_rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                                 "Birth Year","year","tillerid"))
# View(elri_rseedmerge_ss)

elri_rseedmerge_ssf <- merge(elri_rseedmerge_ss, elri_rflw, by = c( "plot", "pos", "tag", "Endo", 
                                                                    "Birth Year","year"), all = TRUE)
# View(elri_rseedmerge_ssf)






# Combining the  original and recruit ELRI repro dataframes ---------
elri_seedmerge_ssf <- elri_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                           "tillerid", "flw","spikelets", "seed")]

elri_rseedmerge_ssf <- elri_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                             "tillerid", "flw","spikelets", "seed")]


ELRIrepro <- elri_seedmerge_ssf %>% 
  rbind(elri_rseedmerge_ssf) %>% 
  mutate(species = "ELRI")

# ELRIrepro <- ELRIrepro[!(is.na(ELRIrepro$flw)),]

# View(ELRIrepro)







# Combining repro measurements across years for the AGPE data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken
### The repro data for AGPE is pretty funky. These are recorded as averages from multiple tillers as opposed to as counts with tiller id's for the other species.
### The 2016 spikelet values are total spikelets for the plants, but all but one of the plants have only one tiller, meaning that for most of the plants the data is essentially spikelets per inflorescence. Currently I left this in, but the plant is plot 120, Pos 17, tag 2397. 
### spikelet data is also recorded as tot spikes for a few years of recruits data.
agpeseed <- AGPE_data %>%
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(seed2007 = NA,  seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seeds_spikelet1", "seeds_spikelet2", 
                       "seeds_spikelet3", "seed2011",  "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = "multitillermean")
agpeseed$year<-  ifelse(agpeseed$variable == "seed2007", 2007, ifelse(agpeseed$variable == "seeds_spikelet1", 2008, ifelse(agpeseed$variable  == "seeds_spikelet2", 2009, ifelse(agpeseed$variable  == "seeds_spikelet3", 2010, ifelse(agpeseed$variable  == "seed2011", 2011, ifelse(agpeseed$variable  == "seed2012", 2012, ifelse(agpeseed$variable  == "seed2013", 2013,ifelse(agpeseed$variable == "seed2014", 2014,ifelse(agpeseed$variable == "seed2015", 2015,ifelse(agpeseed$variable  == "seed2016", 2016, NA))))))))))
agpeseed1 <- agpeseed %>% 
  filter(!is.na(seed))
# View(agpeseed1)

agpespike <- AGPE_data %>% 
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(spike2007 = NA,) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "Spikelets_tiller1", "Spikelets_Infl2", "no_total_spikelets_infl3",
                       "avg_spikelets4", "avg_spikelets5", "avg_spikelets6", 
                       "avg_spikelets7", "spike_infl8", "TotSpikelets9"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = "multitillermean")
agpespike$year<- ifelse(agpespike$variable == "spike2007", 2007, ifelse(agpespike$variable == "Spikelets_tiller1", 2008, ifelse(agpespike$variable  == "Spikelets_Infl2", 2009, ifelse(agpespike$variable  == "no_total_spikelets_infl3", 2010, ifelse(agpespike$variable  == "avg_spikelets4", 2011, ifelse(agpespike$variable  == "avg_spikelets5", 2012, ifelse(agpespike$variable  == "avg_spikelets6", 2013, ifelse(agpespike$variable == "avg_spikelets7", 2014, ifelse(agpespike$variable == "spike_infl8", 2015, ifelse(agpespike$variable  == "TotSpikelets9", 2016, NA))))))))))
agpespike1 <- agpespike %>% 
  filter(!is.na(spikelets))
# View(agpespike1)
 

agpeflw <- AGPE_data %>% 
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FLWTiller3", 
                       "FLWTiller4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
agpeflw$year<- ifelse(agpeflw$variable == "flw2007", 2007, ifelse(agpeflw$variable == "Flwtillers1", 2008, ifelse(agpeflw$variable  == "FlwTillers2", 2009, ifelse(agpeflw$variable  == "FLWTiller3", 2010, ifelse(agpeflw$variable  == "FLWTiller4", 2011, ifelse(agpeflw$variable  == "FlwTillers5", 2012, ifelse(agpeflw$variable  == "FlwTillers6", 2013,ifelse(agpeflw$variable == "FlwTillers7", 2014,ifelse(agpeflw$variable == "FlwTillers8", 2015,ifelse(agpeflw$variable  == "FlwTillers9", 2016, NA))))))))))
# View(agpeflw)

agpe_seedmerge_ss <- full_join(agpespike1, agpeseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                    "Birth Year","year","tillerid"))
# View(agpe_seedmerge_ss)

agpe_seedmerge_ssf <- merge(agpe_seedmerge_ss, agpeflw, by = c("plot", "pos", "tag", "Endo", 
                                                       "Birth Year","year"), all = TRUE)
# View(agpe_seedmerge_ssf)


# Combining repro measurements across years for the AGPE recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
### The spikelets data is also reported in Total spikelets for several years. The first two years there aren't any flwing tillers, and then in 2011 there are only single tillers, so this could essentially be spikelets/infl
### 2013 has lots of spikelet data, but this it is recorded as totals for the plant within an equatio, so I am calculating avgs from it with the number of flw tillers
### 2014 has data recorded as avg spikelets
agpe_rseed <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed")  %>% 
  mutate(tillerid = "multitillermean")
agpe_rseed$year<- ifelse(agpe_rseed$variable == "seed2007", 2007, ifelse(agpe_rseed$variable == "seed2008", 2008, ifelse(agpe_rseed$variable == "seed2009", 2009, ifelse(agpe_rseed$variable == "seed2010", 2010, ifelse(agpe_rseed$variable == "seed2011", 2011, ifelse(agpe_rseed$variable  == "seed2012", 2012, ifelse(agpe_rseed$variable  == "seed2013", 2013, ifelse(agpe_rseed$variable  == "seed2014", 2014, ifelse(agpe_rseed$variable  == "seed2015", 2015, ifelse(agpe_rseed$variable == "seed2016", 2016, NA))))))))))
agpe_rseed1 <- agpe_rseed %>% 
  filter(!is.na(seed))
# View(agpe_rseed1)

agpe_rspike <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA) %>% 
  mutate(avgspike13 = TotSpikelets13/FlwTillers13) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", 
                       "spike2010","TotSpikelets11", "SpikeletsA12", 
                       "avgspike13", "avgspikepertiller14", "spikepertillerA15","spikepertillerB15",
                       "spikepertillerA16","spikepertillerB16"),
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("avg", variable) ~ "multitillermean",
                              grepl("Spikelets", variable) ~ "multitillermean"))
agpe_rspike$year<-  ifelse(agpe_rspike$variable == "spike2007", 2007, ifelse(agpe_rspike$variable == "spike2008", 2008, ifelse(agpe_rspike$variable == "spike2009", 2009, ifelse(agpe_rspike$variable == "spike2010", 2010, ifelse(agpe_rspike$variable == "TotSpikelets11", 2011, ifelse(agpe_rspike$variable  == "SpikeletsA12", 2012, ifelse(agpe_rspike$variable  == "avgspike13", 2013, ifelse(agpe_rspike$variable  == "avgspikepertiller14", 2014, ifelse(agpe_rspike$variable  == "spikepertillerA15", 2015, ifelse(agpe_rspike$variable  == "spikepertillerB15", 2015,ifelse(agpe_rspike$variable == "spikepertillerA16", 2016, ifelse(agpe_rspike$variable == "spikepertillerB16", 2016, NA))))))))))))
agpe_rspike1 <- agpe_rspike %>% 
  filter(!is.na(spikelets))
# View(agpe_rspike1)

agpe_rflw <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(flw2007 = NA, flw2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "FlwTillers09", "FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15", "FlwTillers16"),
       value.name = "flw") 
agpe_rflw$year<- ifelse(agpe_rflw$variable == "flw2007", 2007, ifelse(agpe_rflw$variable == "flw2008", 2008, ifelse(agpe_rflw$variable == "FlwTillers09", 2009, ifelse(agpe_rflw$variable == "FlwTillers10", 2010, ifelse(agpe_rflw$variable == "FlwTillers11", 2011, ifelse(agpe_rflw$variable  == "FlwTillers12", 2012, ifelse(agpe_rflw$variable  == "FlwTillers13", 2013, ifelse(agpe_rflw$variable  == "FlwTillers14", 2014, ifelse(agpe_rflw$variable  == "FlwTillers15", 2015, ifelse(agpe_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(agpe_rflw)

agpe_rseedmerge_ss <- full_join(agpe_rspike1, agpe_rseed1, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year", "tillerid"))
# View(agpe_rseedmerge_ss)

agpe_rseedmerge_ssf <- merge(agpe_rseedmerge_ss, agpe_rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"), all = TRUE)
# View(agpe_rseedmerge_ssf)






# Combining the  original and recruit AGPE repro dataframes ---------
agpe_seedmerge_ssf <- agpe_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                           "tillerid", "flw","spikelets", "seed")]

agpe_rseedmerge_ssf <- agpe_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                             "tillerid", "flw","spikelets", "seed")]


AGPErepro <- agpe_seedmerge_ssf %>% 
  rbind(agpe_rseedmerge_ssf) %>% 
  mutate(species = "AGPE")



# View(AGPErepro)

###### Bind together the different species repro datasets and merge with endo_demog_long
LTREB_repro <- AGPErepro %>% 
  merge(ELRIrepro, by = c("plot", "pos","tag", "Endo", "Birth Year", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(ELVIrepro, by = c("plot", "pos","tag", "Endo", "Birth Year", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(FESUrepro, by = c("plot", "pos","tag", "Endo", "Birth Year", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(LOARrepro, by = c("plot", "pos","tag", "Endo", "Birth Year", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(POALrepro, by = c("plot", "pos","tag", "Endo", "Birth Year", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(POSYrepro, by = c("plot", "pos","tag", "Endo", "Birth Year", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  mutate(flw = as.numeric(flw), seed = as.numeric(seed), spikelets = as.numeric(spikelets))
  
# There are a few instances within the reproductive data where there is a 0 flw measurement with spikelet data. Often these are recorded in the raw data with a note about how seeds were possibly mislabeled.
# There are more instances where there is a flw measurement but spikelet data is not recorded. In some instances, this is likely due to forgetting to collect the data, and there are sometimes in the raw data 
# files where I tried to extract the spikelet count and there were only opaque seed estimates and no specific spikelet counts recorded (although they possibly were at some point to generate the seed estimate). 
# An example of the latter would be POAL plot 3, tag 41, in year 2013, where the note marks that this tiller was skipped.
# I'm going to filter out cases where seeds or spikelet > 0 and flw = 0
LTREB_repro_flw_spike_mismatches <- LTREB_repro %>% 
  filter(flw==0 & spikelets>0 | flw > 0 & is.na(spikelets))

###########################################################################################################################################################################
###### Cleaning up seed and spikelet information and merging with LTREB_data that will be used in the reproductive kernels -----------------
###########################################################################################################################################################################
LTREB_repro1 <- LTREB_repro %>% 
  rename(endo = Endo) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  mutate(`birth` = as.integer(`Birth Year`)) %>% 
  mutate(plot_fixed = as.integer(plot)) %>% 
  mutate(spikelets_fixed = as.numeric(case_when(flw > 0 & !is.na(spikelets) ~ as.character(spikelets),
                                     flw == 0 & spikelets == 0 ~ NA_character_,
                                     flw == 0 & spikelets != 0 ~ NA_character_,
                                     flw == 0 & is.na(spikelets) ~ NA_character_))) %>% 
  mutate(tillerid_fixed = case_when(!is.na(tillerid) & spikelets >= 0 ~ tillerid,
                                    !is.na(tillerid) & seed >= 0 ~ tillerid,
                               is.na(tillerid) & !is.na(spikelets) ~ "A",
                               is.na(tillerid) & is.na(spikelets) ~ NA_character_)) %>% 
  distinct()

# spreading out the spikelet info by tiller to create single row per year per individual.
LTREB_repro2 <- LTREB_repro1 %>%
  select(plot_fixed, pos, tag, endo_01, birth, year, species, flw, spikelets_fixed, tillerid_fixed) %>% 
  spread(key = tillerid_fixed, value = spikelets_fixed) %>% 
  rename(spike_a_t1 = A, spike_b_t1 = B, spike_c_t1 = C, spike_d_t1 = D,  spikelets_AGPE_mean = multitillermean) %>% 
  select(-'<NA>')
  
  
# View(LTREB_repro2)
dim(LTREB_repro2)
table(LTREB_repro2$species, LTREB_repro2$year, !is.na(LTREB_repro2$flw))
table(is.na(LTREB_repro1$flw), LTREB_repro1$flw>0, !is.na(LTREB_repro1$spikelets))

# I am going to try to merge and then add the lagged repro variables at the end
LTREB_repro2_t1 <- LTREB_repro2%>% 
  rename(flw_t1 = flw,
         year_t1 = year)


# merge the reproductive data with LTREB long data file for recent data and for size information for the reproductive model
# This is endodemoglong which is stored in LTREB_data
head(LTREB_data)

LTREB_data1 <- LTREB_data %>% 
  select(-contains("seed"), -plot, -endo)



# now we can merge the two datasets together.
LTREB_repro_combo <- LTREB_data1 %>% 
  left_join(LTREB_repro2_t1,
            by = c("plot_fixed" = "plot_fixed", "pos" = "pos",
                   "id" = "tag", "species" = "species",
                   "endo_01" = "endo_01",
                   "year_t1" = "year_t1")) %>% 
  rename("birth" = "birth.x", "birth_fromrepro" = "birth.y",
         "spike_a_t1_long" = "spike_a_t1.x", "spike_a_t1_fromrepro" = "spike_a_t1.y",
         "spike_b_t1_long" = "spike_b_t1.x", "spike_b_t1_fromrepro" = "spike_b_t1.y",
         "spike_c_t1_long" = "spike_c_t1.x", "spike_c_t1_fromrepro" = "spike_c_t1.y",
         "spike_d_t1_fromrepro" = "spike_d_t1", 
         "flw_t1_long" = "flw_t1.x", "flw_t1_fromrepro" = "flw_t1.y")
# View(LTREB_repro_combo)



# Now select the correct repro data from long or from the raw files into new master columns
LTREB_full_to2018 <- LTREB_repro_combo %>% 
  mutate(FLW_COUNT_T1 = as.integer(case_when(surv_t1 == 0 ~ NA_integer_, !is.na(flw_t1_fromrepro) & !is.na(flw_t1_long) ~ as.integer(flw_t1_fromrepro),
                                       is.na(flw_t1_fromrepro) & !is.na(flw_t1_long) ~ as.integer(flw_t1_long),
                                       !is.na(flw_t1_fromrepro) & is.na(flw_t1_long) ~ as.integer(flw_t1_fromrepro))),   # In this case, where both datasets had data entered, spot checking showed that they have they were identical
         FLW_STAT_T1 = as.integer(case_when(FLW_COUNT_T1 > 0 & surv_t1 == 1 ~ 1, FLW_COUNT_T1 == 0 & surv_t1 == 1 ~ 0, is.na(FLW_COUNT_T1) & surv_t1 == 1 ~ 0)),
         SPIKE_A_T1 =  case_when(is.na(spike_a_t1_fromrepro) & !is.na(spike_a_t1_long) ~ as.numeric(spike_a_t1_long),
                                 !is.na(spike_a_t1_fromrepro) & is.na(spike_a_t1_long) ~ as.numeric(spike_a_t1_fromrepro),                                                                          
                                 !is.na(spike_a_t1_fromrepro) & !is.na(spike_a_t1_long) ~ as.numeric(spike_a_t1_fromrepro)),
         SPIKE_B_T1 =  case_when(is.na(spike_b_t1_fromrepro) & !is.na(spike_b_t1_long) ~ as.numeric(spike_b_t1_long),
                                 !is.na(spike_b_t1_fromrepro) & is.na(spike_b_t1_long) ~ as.numeric(spike_b_t1_fromrepro),                                                                          
                                 !is.na(spike_b_t1_fromrepro) & !is.na(spike_b_t1_long) ~ as.numeric(spike_b_t1_fromrepro)),
         SPIKE_C_T1 =  case_when(is.na(spike_c_t1_fromrepro) & !is.na(spike_c_t1_long) ~ as.numeric(spike_c_t1_long),
                                 !is.na(spike_c_t1_fromrepro) & is.na(spike_c_t1_long) ~ as.numeric(spike_c_t1_fromrepro),                                                                          
                                 !is.na(spike_c_t1_fromrepro) & !is.na(spike_c_t1_long) ~ as.numeric(spike_c_t1_fromrepro)),
         SPIKE_D_T1 = spike_d_t1_fromrepro)
         

# Then we will add lagged variables to have the measurements in time t
LTREB_full_to2018_lag <- LTREB_full_to2018 %>% 
  group_by(id) %>% 
  mutate(FLW_COUNT_T = dplyr::lag(FLW_COUNT_T1, n = 1, default = NA),
         FLW_STAT_T = dplyr::lag(FLW_STAT_T1, n = 1, default = NA),
         SPIKE_A_T = dplyr::lag(SPIKE_A_T1, n = 1, default = NA),
         SPIKE_B_T = dplyr::lag(SPIKE_B_T1, n = 1, default = NA),
         SPIKE_C_T = dplyr::lag(SPIKE_C_T1, n = 1, default = NA),
         SPIKE_D_T = dplyr::lag(SPIKE_D_T1, n = 1, default = NA)) %>% 
  select(plot_fixed, pos, id, species, species_index, 
         endo_01, endo_index, origin_01, birth,
         year_t1, year_t1_index,
         surv_t1, size_t1, logsize_t1,
         FLW_COUNT_T1, FLW_STAT_T1,
         SPIKE_A_T1, SPIKE_B_T1, SPIKE_C_T1, SPIKE_D_T1,
         year_t, year_t_index, size_t, logsize_t, 
         FLW_COUNT_T, FLW_STAT_T,
         SPIKE_A_T, SPIKE_B_T, SPIKE_C_T, SPIKE_D_T)

##############################################################################
####### Here we will merge in 2019 data ------------------------------
##############################################################################

# ELRI_2019_data <- read_xlsx("~/Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019_July2019FieldData_NumbersConversion.xlsx", sheet = "ELRI")
# ELVI_2019_data <- read_xlsx("~/Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019_July2019FieldData_NumbersConversion.xlsx", sheet = "ELVI")
# FESU_2019_data <- read_xlsx("~/Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019_July2019FieldData_NumbersConversion.xlsx", sheet = "FESU")
# LOAR_2019_data <- read_xlsx("~/Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019_July2019FieldData_NumbersConversion.xlsx", sheet = "LOAR")
# POAL_2019_data <- read_xlsx("~/Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019_July2019FieldData_NumbersConversion.xlsx", sheet = "POAL") %>% 
#   mutate(id = paste(plot, "_", pos))
# POSY_2019_data <- read_xlsx("~/Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019_July2019FieldData_NumbersConversion.xlsx", sheet = "POSY") %>% 
#   mutate(id = paste(plot, "_", pos))
# 
# # Now we can merge all the different species together.
# LTREB_2019_data <- ELRI_2019_data %>% 
#   merge(ELVI_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
#   merge(FESU_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
#   merge(LOAR_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
#   merge(POAL_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
#   merge(POSY_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE)
#   
# # We need to do a little bit of cleaning up for some of the missing tags and changing variable names.
# LTREB_2019_cleaned <- LTREB_2019_data %>% 
#   rename(year_t1 = observation_year, surv_t1 = survival, 
#          size_t1 = size_tillers, FLW_COUNT_T1 = flowering_tillers, 
#          SPIKE_A_T1 = spikelets_A, SPIKE_B_T1 = spikelets_B, SPIKE_C_T1 = spikelets_C,
#          dist_a = distance_A, dist_b = distance_B)
#   mutate(plot_fixed = as.integer(plot))
#   
#   
#   # I'm gonna talk to Tom about this stuff and see if the excel sheet can be cleaned up more.
#   
#   LTREB_data <- LTREB_endodemog %>% 
#     mutate(size_t = na_if(size_t, 0)) %>% 
#     mutate(size_t1 = na_if(size_t1, 0)) %>%  
#     mutate(size_t, logsize_t = log(size_t)) %>% 
#     mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
#     mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
#     mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
#                                           endo == "1"| endo =="plus" ~ 1))) %>% 
#     mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
#     mutate(species = case_when(species == "ELVI" & plot == 101 ~ "ELRI", species == "ELVI" & plot != 101 ~ "ELVI",  # This is for the compiled data where a ELRI data point in plot 101, tag 2004 is labelled as ELVI
#                                species == "ELRI" ~ "ELRI",
#                                species == "FESU" ~ "FESU",
#                                species == "AGPE" ~ "AGPE",
#                                species == "POAL" ~ "POAL",
#                                species == "POSY" ~ "POSY",
#                                species == "LOAR" ~ "LOAR"))  %>%    
#     mutate(species_index = as.integer(recode_factor(species,                   
#                                                     "AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
#                                                     "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
#                                                     "POSY" = 7))) %>% 
#     mutate(year_t_index = as.integer(recode(year_t, 
#                                             '2007' = 1, '2008' = 2, '2009' = 3, 
#                                             '2010' = 4, '2011' = 5, '2012' = 6, 
#                                             '2013' = 7, '2014' = 8, '2015' = 9, 
#                                             '2016' = 10, '2017' = 11))) %>%             
#     mutate(year_t1_index = as.integer(recode(year_t1, 
#                                              '2008' = 2, '2009' = 3, '2010' = 4, 
#                                              '2011' = 5, '2012' = 6, '2013' = 7, 
#                                              '2014' = 8, '2015' = 9, '2016' = 10, 
#                                              '2017' = 11, '2018' = 12))) %>%               
#     mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, 
#                                             origin == "R" ~ 1, 
#                                             origin != "R" | origin != "O" ~ 1))) %>%
#     mutate(plot_fixed = as.integer(case_when(species == "LOAR" & id == "39_1B" ~ "39", # This is for a copy error with LOAR individual 39_1B, which assigned it to plots 41-44
#                                              plot != "R" ~ as.character(plot), 
#                                              plot == "R" ~ as.character(origin)))) %>% 
#     mutate(surv_t1 = as.integer(case_when(surv_t1 == 1 ~ 1,
#                                           surv_t1 == 0 ~ 0,
#                                           is.na(surv_t1) & birth == year_t1 ~ 1))) %>% 
#     filter(duplicated(.) == FALSE)
#   # dim(LTREB_data)
#   
#   
#   
# 
# # We still need 2019 AGPE data
# LTREB_full_to2019 <- LTREB_full_to2018_lag %>% 
#   full_join()
# 

##############################################################################
####### Merging in the endophyte checks ------------------------------
##############################################################################

LTREB_endo_check <- read_csv(file = "~/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_status.csv") %>% 
  select(-recno,-check, -X11) %>% 
  rename("origin_from_check" = "origin", "endo_status_from_check" = "status", "plot_endo_for_check" = "endo") %>% 
  mutate(origin_01 = as.integer(case_when(origin_from_check == "O" ~ 0, 
                                          origin_from_check == "R" ~ 1))) %>% 
  mutate(plot_endo_for_check = as.integer(recode(plot_endo_for_check, "plus" = 1, "minus" = 0))) %>% 
  mutate(endo_mismatch = plot_endo_for_check - endo_status_from_check) # 0 = no change, >0 = loss of endophyte from positive plot, <0 = gain of endophyte in negative plot
  # Metadata from Jenn
#   recno:	unique record number
#   species:	four letter species code
#   origin:	"O" = original plant from greenhouse, "R" = recruit
#   plot:	plot number
#   pos:	for "O" = position of plant in plot, for "R" = recruit number on tag
#   id:	unique identifier for each plant, "O" = single unique number, "R" = concatenation of plot number and tag number separated by "_"
#   status:	0= no endophyte found via microscopy on leaf peel at 200X, 1= endophyte detected
#   date_status:	year leaf peel was taken
#   endo:	plot level endophyte status: "minus" = no endophyte in original planting, "plus" = endophyte present in original planting
#   check:	"x" indicates the incorrect status was detected given original plot level treatment

# There are two plants that are checked but are not present in the endo_demog_long dataset
setdiff(LTREB_endo_check$id,LTREB_full_to2018_lag$id)

LTREB_full_2 <- LTREB_full_to2018_lag %>% 
  left_join(LTREB_endo_check, by = c("species" = "species", "plot_fixed" = "plot", "pos" = "pos", "origin_01" = "origin_01", "id" = "id"))

# here are some summaries of the amount of changes in endophyte status
LTREB_status_changes <- LTREB_full_2 %>% 
  group_by(species, plot_fixed) %>% 
  summarize(same = sum(endo_mismatch == 0, na.rm = TRUE),
            lose_endo = sum(endo_mismatch > 0, na.rm = TRUE),
            gain_endo = sum(endo_mismatch <0, na.rm = TRUE)) %>% 
  mutate(percent_gain = (gain_endo/same)*100, percent_lose = (lose_endo/same)*100)
  

##############################################################################
####### Merging in the location data ------------------------------
##############################################################################

LTREB_distances <- read_csv(file = "~/Dropbox/EndodemogData/Fulldataplusmetadata/endo_distance_tubeid.csv",
                            col_types = cols( species = col_character(),
                            origin = col_character(),
                            plot = col_integer(),
                            pos = col_character(),
                            id = col_character(),
                            dist_a = col_double(),
                            dist_b = col_double(),
                            tubeid = col_character(),
                            notes = col_character(),
                            date_dist = col_character())) %>% 
  select(species, origin, plot, pos, id, dist_a, dist_b, date_dist) %>% 
  rename("origin_from_distance" = "origin") %>% 
  mutate(origin_01 = as.integer(case_when(origin_from_distance == "O" ~ 0, 
                                          origin_from_distance == "R" ~ 1))) %>% 
  filter(!is.na(dist_a), !is.na(dist_b)) %>% 
  mutate(duplicate = duplicated(id)) %>% #There are three LOAR id's that have two measurements, one from may 2018 and one from sept 2018: 40_F5, 33_4B, 33_12
  filter(!(species == "LOAR" & date_dist == "may_18" & id %in% c("40_F5", "33_4B", "33_12"))) # I am selecting the september measurements for these id's. The distances are different but similar

# Here are the plant id's that are in the distance file but not the long file
setdiff(LTREB_distances$id, LTREB_full_2$id)


LTREB_full <- LTREB_full_2 %>% 
  left_join(LTREB_distances, by = c("species" = "species","pos" = "pos", "plot_fixed" = "plot", "origin_01" = "origin_01", "id" = "id")) %>% 
  select(-duplicate, -origin_from_check, -origin_from_distance, -date_status, -date_dist) # I'm removing some of the extrneous variable. We also have distance data in the new field data that needs to be merged in.


##############################################################################
####### Preparing datalists for Survival Kernel ------------------------------
##############################################################################

# NA's in survival come from mostly 2017 recruits.
LTREB_data_forsurv <- LTREB_full %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))

dim(LTREB_data_forsurv)

# Creating individual species data lists to be passed to the model

# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"111"="1", "112"="2", "113"="3", "114"="4", "115"="5", "116"="6", "117"="7", "118"="8","119"="9", "120"="10"))))
ELRI_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"101"="1", "102"="2", "103"="3", "104"="4", "105"="5", "106"="6", "107"="7", "108"="8","109"="9", "110"="10"))))
ELVI_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"91"="1", "92"="2", "93"="3", "94"="4", "95"="5", "96"="6", "97"="7", "98"="8","99"="9", "100"="10"))))
FESU_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"121"="1", "122"="2", "123"="3", "124"="4", "125"="5", "126"="6", "127"="7", "128"="8","129"="9", "130"="10"))))
LOAR_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"31"="1", "32"="2", "33"="3", "34"="4", "35"="5", "36"="6", "37"="7", "38"="8","39"="9", "40"="10"))))
POAL_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"3"="1", "4"="2", "8"="3", "9"="4", "10"="5", "11"="6", "15"="7", "16"="8","17"="9", "19"="10","151"="11","152"="12","153"="13","154"="14","155"="15","156"="16","157"="17","158"="18")))) 
POSY_surv_data <- LTREB_data_forsurv %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"1"="1", "2"="2", "5"="3", "6"="4", "7"="5", "12"="6", "13"="7", "14"="8","18"="9", "20"="10","141"="11","142"="12","143"="13","144"="14","145"="15","146"="16","147"="17","148"="18", "149"="19", "150"="20"))))

# Create model matrices for each species year*endo and plot random effects
AGPE_yearendo_Xs_s <- AGPE_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

AGPE_plot_Xs_s <- AGPE_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELRI_yearendo_Xs_s <- ELRI_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

ELRI_plot_Xs_s <- ELRI_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELVI_yearendo_Xs_s <- ELVI_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

ELVI_plot_Xs_s <- ELVI_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


FESU_yearendo_Xs_s <- FESU_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

FESU_plot_Xs_s <- FESU_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

LOAR_yearendo_Xs_s <- LOAR_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

LOAR_plot_Xs_s <- LOAR_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

POAL_yearendo_Xs_s <- POAL_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

POAL_plot_Xs_s <- POAL_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)
  

POSY_yearendo_Xs_s <- POSY_surv_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

POSY_plot_Xs_s <- POSY_surv_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


# Create data lists to be used for the Stan model
AGPE_surv_data_list <- list(surv_t1 = AGPE_surv_data$surv_t1,
                            yearendo_Xs = AGPE_yearendo_Xs_s,
                            plot_Xs = AGPE_plot_Xs_s,
                            logsize_t = AGPE_surv_data$logsize_t,
                            origin_01 = AGPE_surv_data$origin_01,
                            endo_01 = AGPE_surv_data$endo_01,
                            endo_index = AGPE_surv_data$endo_index,
                            year_t = AGPE_surv_data$year_t_index,
                            plot = AGPE_surv_data$plot_index,
                            N = nrow(AGPE_surv_data),
                            K = 5L,
                            nYear = length(unique(AGPE_surv_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(AGPE_surv_data$endo_01)))
str(AGPE_surv_data_list)

ELRI_surv_data_list <- list(surv_t1 = ELRI_surv_data$surv_t1,
                            yearendo_Xs = ELRI_yearendo_Xs_s,
                            plot_Xs = ELRI_plot_Xs_s,
                            logsize_t = ELRI_surv_data$logsize_t,
                            origin_01 = ELRI_surv_data$origin_01,
                            endo_01 = ELRI_surv_data$endo_01,
                            endo_index = ELRI_surv_data$endo_index,
                            year_t = ELRI_surv_data$year_t_index,
                            plot = ELRI_surv_data$plot_index,
                            N = nrow(ELRI_surv_data),
                            K = 5L,
                            nYear = length(unique(ELRI_surv_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(ELRI_surv_data$endo_01)))
str(ELRI_surv_data_list)

ELVI_surv_data_list <- list(surv_t1 = ELVI_surv_data$surv_t1,
                            yearendo_Xs = ELVI_yearendo_Xs_s,
                            plot_Xs = ELVI_plot_Xs_s,
                            logsize_t = ELVI_surv_data$logsize_t,
                            origin_01 = ELVI_surv_data$origin_01,
                            endo_01 = ELVI_surv_data$endo_01,
                            endo_index = ELVI_surv_data$endo_index,
                            year_t = ELVI_surv_data$year_t_index,
                            plot = ELVI_surv_data$plot_index,
                            N = nrow(ELVI_surv_data),
                            K = 5L,
                            nYear = length(unique(ELVI_surv_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(ELVI_surv_data$endo_01)))
str(ELVI_surv_data_list)

FESU_surv_data_list <- list(surv_t1 = FESU_surv_data$surv_t1,
                            yearendo_Xs = FESU_yearendo_Xs_s,
                            plot_Xs = FESU_plot_Xs_s,
                            logsize_t = FESU_surv_data$logsize_t,
                            origin_01 = FESU_surv_data$origin_01,
                            endo_01 = FESU_surv_data$endo_01,
                            endo_index = FESU_surv_data$endo_index,
                            year_t = FESU_surv_data$year_t_index,
                            plot = FESU_surv_data$plot_index,
                            N = nrow(FESU_surv_data),
                            K = 5L,
                            nYear = length(unique(FESU_surv_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(FESU_surv_data$endo_01)))
str(FESU_surv_data_list)

LOAR_surv_data_list <- list(surv_t1 = LOAR_surv_data$surv_t1,
                            yearendo_Xs = LOAR_yearendo_Xs_s,
                            plot_Xs = LOAR_plot_Xs_s,
                            logsize_t = LOAR_surv_data$logsize_t,
                            origin_01 = LOAR_surv_data$origin_01,
                            endo_01 = LOAR_surv_data$endo_01,
                            endo_index = LOAR_surv_data$endo_index,
                            year_t = LOAR_surv_data$year_t_index,
                            plot = LOAR_surv_data$plot_index,
                            N = nrow(LOAR_surv_data),
                            K = 5L,
                            nYear = length(unique(LOAR_surv_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(LOAR_surv_data$endo_01)))
str(LOAR_surv_data_list)

POAL_surv_data_list <- list(surv_t1 = POAL_surv_data$surv_t1,
                            yearendo_Xs = POAL_yearendo_Xs_s,
                            plot_Xs = POAL_plot_Xs_s,
                            logsize_t = POAL_surv_data$logsize_t,
                            origin_01 = POAL_surv_data$origin_01,
                            endo_01 = POAL_surv_data$endo_01,
                            endo_index = POAL_surv_data$endo_index,
                            year_t = POAL_surv_data$year_t_index,
                            plot = POAL_surv_data$plot_index,
                            N = nrow(POAL_surv_data),
                            K = 5L,
                            nYear = length(unique(POAL_surv_data$year_t_index)),
                            nPlot = 18L,
                            nEndo =   length(unique(POAL_surv_data$endo_01)))
str(POAL_surv_data_list)

POSY_surv_data_list <- list(surv_t1 = POSY_surv_data$surv_t1,
                            yearendo_Xs = POSY_yearendo_Xs_s,
                            plot_Xs = POSY_plot_Xs_s,
                            logsize_t = POSY_surv_data$logsize_t,
                            origin_01 = POSY_surv_data$origin_01,
                            endo_01 = POSY_surv_data$endo_01,
                            endo_index = POSY_surv_data$endo_index,
                            year_t = POSY_surv_data$year_t_index,
                            plot = POSY_surv_data$plot_index,
                            N = nrow(POSY_surv_data),
                            K = 5L,
                            nYear = length(unique(POSY_surv_data$year_t_index)),
                            nPlot = 20L,
                            nEndo =   length(unique(POSY_surv_data$endo_01)))
str(POSY_surv_data_list)



##############################################################################
####### Preparing datalists for Growth Kernel ------------------------------
##############################################################################

LTREB_data_forgrow <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(size_t1)) %>% 
  filter(!is.na(endo_01)) %>% 
  mutate(size_t1 = as.integer(size_t1))

dim(LTREB_data_forgrow)

# Creating individual species dataframes and updating plots for indexing

AGPE_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"111"="1", "112"="2", "113"="3", "114"="4", "115"="5", "116"="6", "117"="7", "118"="8","119"="9", "120"="10"))))
ELRI_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"101"="1", "102"="2", "103"="3", "104"="4", "105"="5", "106"="6", "107"="7", "108"="8","109"="9", "110"="10"))))
ELVI_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"91"="1", "92"="2", "93"="3", "94"="4", "95"="5", "96"="6", "97"="7", "98"="8","99"="9", "100"="10"))))
FESU_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"121"="1", "122"="2", "123"="3", "124"="4", "125"="5", "126"="6", "127"="7", "128"="8","129"="9", "130"="10"))))
LOAR_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"31"="1", "32"="2", "33"="3", "34"="4", "35"="5", "36"="6", "37"="7", "38"="8","39"="9", "40"="10"))))
POAL_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"3"="1", "4"="2", "8"="3", "9"="4", "10"="5", "11"="6", "15"="7", "16"="8","17"="9", "19"="10","151"="11","152"="12","153"="13","154"="14","155"="15","156"="16","157"="17","158"="18")))) 
POSY_grow_data <- LTREB_data_forgrow %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"1"="1", "2"="2", "5"="3", "6"="4", "7"="5", "12"="6", "13"="7", "14"="8","18"="9", "20"="10","141"="11","142"="12","143"="13","144"="14","145"="15","146"="16","147"="17","148"="18", "149"="19", "150"="20"))))






# Create model matrices for each species year*endo and plot random effects
AGPE_yearendo_Xs_g <- AGPE_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

AGPE_plot_Xs_g <- AGPE_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELRI_yearendo_Xs_g <- ELRI_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

ELRI_plot_Xs_g <- ELRI_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELVI_yearendo_Xs_g <- ELVI_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

ELVI_plot_Xs_g <- ELVI_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


FESU_yearendo_Xs_g <- FESU_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

FESU_plot_Xs_g <- FESU_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

LOAR_yearendo_Xs_g <- LOAR_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

LOAR_plot_Xs_g <- LOAR_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

POAL_yearendo_Xs_g <- POAL_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

POAL_plot_Xs_g <- POAL_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


POSY_yearendo_Xs_g <- POSY_grow_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)

POSY_plot_Xs_g <- POSY_grow_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

# Create data lists to be used for the Stan model
AGPE_grow_data_list <- list(size_t1 = AGPE_grow_data$size_t1,
                            yearendo_Xs = AGPE_yearendo_Xs_g,
                            plot_Xs = AGPE_plot_Xs_g,
                            logsize_t = AGPE_grow_data$logsize_t,
                            origin_01 = AGPE_grow_data$origin_01,
                            endo_01 = AGPE_grow_data$endo_01,
                            endo_index = AGPE_grow_data$endo_index,
                            year_t = AGPE_grow_data$year_t_index,
                            plot = AGPE_grow_data$plot_index,
                            N = nrow(AGPE_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(AGPE_grow_data$size_t1)),
                            nYear = length(unique(AGPE_grow_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(AGPE_grow_data$endo_01)))
str(AGPE_grow_data_list)

ELRI_grow_data_list <- list(size_t1 = ELRI_grow_data$size_t1,
                            yearendo_Xs = ELRI_yearendo_Xs_g,
                            plot_Xs = ELRI_plot_Xs_g,
                            logsize_t = ELRI_grow_data$logsize_t,
                            origin_01 = ELRI_grow_data$origin_01,
                            endo_01 = ELRI_grow_data$endo_01,
                            endo_index = ELRI_grow_data$endo_index,
                            year_t = ELRI_grow_data$year_t_index,
                            plot = ELRI_grow_data$plot_index,
                            N = nrow(ELRI_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(ELRI_grow_data$size_t1)),
                            nYear = length(unique(ELRI_grow_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(ELRI_grow_data$endo_01)))
str(ELRI_grow_data_list)

ELVI_grow_data_list <- list(size_t1 = ELVI_grow_data$size_t1,
                            yearendo_Xs = ELVI_yearendo_Xs_g,
                            plot_Xs = ELVI_plot_Xs_g,
                            logsize_t = ELVI_grow_data$logsize_t,
                            origin_01 = ELVI_grow_data$origin_01,
                            endo_01 = ELVI_grow_data$endo_01,
                            endo_index = ELVI_grow_data$endo_index,
                            year_t = ELVI_grow_data$year_t_index,
                            plot = ELVI_grow_data$plot_index,
                            N = nrow(ELVI_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(ELVI_grow_data$size_t1)),
                            nYear = length(unique(ELVI_grow_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(ELVI_grow_data$endo_01)))
str(ELVI_grow_data_list)

FESU_grow_data_list <- list(size_t1 = FESU_grow_data$size_t1,
                            yearendo_Xs = FESU_yearendo_Xs_g,
                            plot_Xs = FESU_plot_Xs_g,
                            logsize_t = FESU_grow_data$logsize_t,
                            origin_01 = FESU_grow_data$origin_01,
                            endo_01 = FESU_grow_data$endo_01,
                            endo_index = FESU_grow_data$endo_index,
                            year_t = FESU_grow_data$year_t_index,
                            plot = FESU_grow_data$plot_index,
                            N = nrow(FESU_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(FESU_grow_data$size_t1)),
                            nYear = length(unique(FESU_grow_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(FESU_grow_data$endo_01)))
str(FESU_grow_data_list)

LOAR_grow_data_list <- list(size_t1 = LOAR_grow_data$size_t1,
                            yearendo_Xs = LOAR_yearendo_Xs_g,
                            plot_Xs = LOAR_plot_Xs_g,
                            logsize_t = LOAR_grow_data$logsize_t,
                            origin_01 = LOAR_grow_data$origin_01,
                            endo_01 = LOAR_grow_data$endo_01,
                            endo_index = LOAR_grow_data$endo_index,
                            year_t = LOAR_grow_data$year_t_index,
                            plot = LOAR_grow_data$plot_index,
                            N = nrow(LOAR_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(LOAR_grow_data$size_t1)),
                            nYear = length(unique(LOAR_grow_data$year_t_index)),
                            nPlot = 10L,
                            nEndo =   length(unique(LOAR_grow_data$endo_01)))
str(LOAR_grow_data_list)

POAL_grow_data_list <- list(size_t1 = POAL_grow_data$size_t1,
                            yearendo_Xs = POAL_yearendo_Xs_g,
                            plot_Xs = POAL_plot_Xs_g,
                            logsize_t = POAL_grow_data$logsize_t,
                            origin_01 = POAL_grow_data$origin_01,
                            endo_01 = POAL_grow_data$endo_01,
                            endo_index = POAL_grow_data$endo_index,
                            year_t = POAL_grow_data$year_t_index,
                            plot = POAL_grow_data$plot_index,
                            N = nrow(POAL_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(POAL_grow_data$size_t1)),
                            nYear = length(unique(POAL_grow_data$year_t_index)),
                            nPlot = 18L,
                            nEndo =   length(unique(POAL_grow_data$endo_01)))
str(POAL_grow_data_list)

POSY_grow_data_list <- list(size_t1 = POSY_grow_data$size_t1,
                            yearendo_Xs = POSY_yearendo_Xs_g,
                            plot_Xs = POSY_plot_Xs_g,
                            logsize_t = POSY_grow_data$logsize_t,
                            origin_01 = POSY_grow_data$origin_01,
                            endo_01 = POSY_grow_data$endo_01,
                            endo_index = POSY_grow_data$endo_index,
                            year_t = POSY_grow_data$year_t_index,
                            plot = POSY_grow_data$plot_index,
                            N = nrow(POSY_grow_data),
                            K = 5L,
                            lowerlimit = as.integer(min(POSY_grow_data$size_t1)),
                            nYear = length(unique(POSY_grow_data$year_t_index)),
                            nPlot = 20L,
                            nEndo =   length(unique(POSY_grow_data$endo_01)))
str(POSY_grow_data_list)


##############################################################################
####### Preparing datalists for Flowering Kernel ------------------------------
##############################################################################

## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data_forflw <- LTREB_full %>% 
  filter(!is.na(FLW_STAT_T)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01))


dim(LTREB_data_forflw)

# Creating individual species data lists to be passed to the model
# Split up the main dataframe by species and recode plot to be used as an index for each species

AGPE_flw_data <- LTREB_data_forflw %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"111"="1", "112"="2", "113"="3", "114"="4", "115"="5", "116"="6", "117"="7", "118"="8","119"="9", "120"="10"))))
ELRI_flw_data <- LTREB_data_forflw %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"101"="1", "102"="2", "103"="3", "104"="4", "105"="5", "106"="6", "107"="7", "108"="8","109"="9", "110"="10"))))
ELVI_flw_data <- LTREB_data_forflw %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"91"="1", "92"="2", "93"="3", "94"="4", "95"="5", "96"="6", "97"="7", "98"="8","99"="9", "100"="10"))))
FESU_flw_data <- LTREB_data_forflw %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"121"="1", "122"="2", "123"="3", "124"="4", "125"="5", "126"="6", "127"="7", "128"="8","129"="9", "130"="10"))))
LOAR_flw_data <- LTREB_data_forflw %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"31"="1", "32"="2", "33"="3", "34"="4", "35"="5", "36"="6", "37"="7", "38"="8","39"="9", "40"="10"))))
POAL_flw_data <- LTREB_data_forflw %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"3"="1", "4"="2", "8"="3", "9"="4", "10"="5", "11"="6", "15"="7", "16"="8","17"="9", "19"="10","151"="11","152"="12","153"="13","154"="14","155"="15","156"="16","157"="17","158"="18"))))
POSY_flw_data <- LTREB_data_forflw %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"1"="1", "2"="2", "5"="3", "6"="4", "7"="5", "12"="6", "13"="7", "14"="8","18"="9", "20"="10","141"="11","142"="12","143"="13","144"="14","145"="15","146"="16","147"="17","148"="18", "149"="19", "150"="20"))))

# Create model matrices for each species year*endo and plot random effects
AGPE_yearendo_Xs_flw <- AGPE_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

AGPE_plot_Xs_flw <- AGPE_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELRI_yearendo_Xs_flw <- ELRI_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)


ELRI_plot_Xs_flw <- ELRI_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELVI_yearendo_Xs_flw <- ELVI_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)%>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

ELVI_plot_Xs_flw <- ELVI_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


FESU_yearendo_Xs_flw <- FESU_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

FESU_plot_Xs_flw <- FESU_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

LOAR_yearendo_Xs_flw <- LOAR_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

LOAR_plot_Xs_flw <- LOAR_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

POAL_yearendo_Xs_flw <- POAL_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

POAL_plot_Xs_flw <- POAL_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


POSY_yearendo_Xs_flw <- POSY_flw_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

POSY_plot_Xs_flw <- POSY_flw_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


# Create data lists to be used for the Stan model
AGPE_flw_data_list <- list(flw_t = AGPE_flw_data$FLW_STAT_T,
                           yearendo_Xs = AGPE_yearendo_Xs_flw,
                           plot_Xs = AGPE_plot_Xs_flw,
                           logsize_t = AGPE_flw_data$logsize_t,
                           origin_01 = AGPE_flw_data$origin_01,
                           endo_01 = AGPE_flw_data$endo_01,
                           endo_index = AGPE_flw_data$endo_index,
                           year_t = AGPE_flw_data$year_t_index,
                           plot = AGPE_flw_data$plot_index,
                           N = nrow(AGPE_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 10L,
                           nEndo =   length(unique(AGPE_flw_data$endo_01)))
str(AGPE_flw_data_list)

ELRI_flw_data_list <- list(flw_t = ELRI_flw_data$FLW_STAT_T,
                           yearendo_Xs = ELRI_yearendo_Xs_flw,
                           plot_Xs = ELRI_plot_Xs_flw,
                           ogsize_t = ELRI_flw_data$logsize_t,
                           origin_01 = ELRI_flw_data$origin_01,
                           endo_01 = ELRI_flw_data$endo_01,
                           endo_index = ELRI_flw_data$endo_index,
                           year_t = ELRI_flw_data$year_t_index,
                           plot = ELRI_flw_data$plot_index,
                           N = nrow(ELRI_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 10L,
                           nEndo =   length(unique(ELRI_flw_data$endo_01)))
str(ELRI_flw_data_list)

ELVI_flw_data_list <- list(flw_t = ELVI_flw_data$FLW_STAT_T,
                           yearendo_Xs = ELVI_yearendo_Xs_flw,
                           plot_Xs = ELVI_plot_Xs_flw,
                           logsize_t = ELVI_flw_data$logsize_t,
                           origin_01 = ELVI_flw_data$origin_01,
                           endo_01 = ELVI_flw_data$endo_01,
                           endo_index = ELVI_flw_data$endo_index,
                           year_t = ELVI_flw_data$year_t_index,
                           plot = ELVI_flw_data$plot_index,
                           N = nrow(ELVI_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 10L,
                           nEndo =   length(unique(ELVI_flw_data$endo_01)))
str(ELVI_flw_data_list)

FESU_flw_data_list <- list(flw_t = FESU_flw_data$FLW_STAT_T,
                           yearendo_Xs = FESU_yearendo_Xs_flw,
                           plot_Xs = FESU_plot_Xs_flw,
                           logsize_t = FESU_flw_data$logsize_t,
                           origin_01 = FESU_flw_data$origin_01,
                           endo_01 = FESU_flw_data$endo_01,
                           endo_index = FESU_flw_data$endo_index,
                           year_t = FESU_flw_data$year_t_index,
                           plot = FESU_flw_data$plot_index,
                           N = nrow(FESU_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 10L,
                           nEndo =   length(unique(FESU_flw_data$endo_01)))
str(FESU_flw_data_list)

LOAR_flw_data_list <- list(flw_t = LOAR_flw_data$FLW_STAT_T,
                           yearendo_Xs = LOAR_yearendo_Xs_flw,
                           plot_Xs = LOAR_plot_Xs_flw,
                           logsize_t = LOAR_flw_data$logsize_t,
                           origin_01 = LOAR_flw_data$origin_01,
                           endo_01 = LOAR_flw_data$endo_01,
                           endo_index = LOAR_flw_data$endo_index,
                           year_t = LOAR_flw_data$year_t_index,
                           plot = LOAR_flw_data$plot_index,
                           N = nrow(LOAR_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 10L,
                           nEndo =   length(unique(LOAR_flw_data$endo_01)))
str(LOAR_flw_data_list)

POAL_flw_data_list <- list(flw_t = POAL_flw_data$FLW_STAT_T,
                           yearendo_Xs = POAL_yearendo_Xs_flw,
                           plot_Xs = POAL_plot_Xs_flw,
                           logsize_t = POAL_flw_data$logsize_t,
                           origin_01 = POAL_flw_data$origin_01,
                           endo_01 = POAL_flw_data$endo_01,
                           endo_index = POAL_flw_data$endo_index,
                           year_t = POAL_flw_data$year_t_index,
                           plot = POAL_flw_data$plot_index,
                           N = nrow(POAL_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 18L,
                           nEndo =   length(unique(POAL_flw_data$endo_01)))
str(POAL_flw_data_list)

POSY_flw_data_list <- list(flw_t = POSY_flw_data$FLW_STAT_T,
                           yearendo_Xs = POSY_yearendo_Xs_flw,
                           plot_Xs = POSY_plot_Xs_flw,
                           logsize_t = POSY_flw_data$logsize_t,
                           origin_01 = POSY_flw_data$origin_01,
                           endo_01 = POSY_flw_data$endo_01,
                           endo_index = POSY_flw_data$endo_index,
                           year_t = POSY_flw_data$year_t_index,
                           plot = POSY_flw_data$plot_index,
                           N = nrow(POSY_flw_data),
                           K = 5L,
                           nYear = 11L,
                           nPlot = 20L,
                           nEndo =   length(unique(POSY_flw_data$endo_01)))
str(POSY_flw_data_list)




##############################################################################
####### Preparing datalists for Fertility Kernel (# of flowering tillers) ------------------------------
##############################################################################


# NA's in survival come from mostly 2017 recruits.
LTREB_data_forfert <- LTREB_full %>% 
  filter(!is.na(FLW_COUNT_T)) %>% 
  filter(FLW_COUNT_T > 0) %>% 
  filter(!is.na(logsize_t))

dim(LTREB_data_forfert)

# Creating individual species data lists to be passed to the model
# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_fert_data <- LTREB_data_forfert %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"111"="1", "112"="2", "113"="3", "114"="4", "115"="5", "116"="6", "117"="7", "118"="8","119"="9", "120"="10"))))
ELRI_fert_data <- LTREB_data_forfert %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"101"="1", "102"="2", "103"="3", "104"="4", "105"="5", "106"="6", "107"="7", "108"="8","109"="9", "110"="10"))))
ELVI_fert_data <- LTREB_data_forfert %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"91"="1", "92"="2", "93"="3", "94"="4", "95"="5", "96"="6", "97"="7", "98"="8","99"="9", "100"="10"))))
FESU_fert_data <- LTREB_data_forfert %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"121"="1", "122"="2", "123"="3", "124"="4", "125"="5", "126"="6", "127"="7", "128"="8","129"="9", "130"="10"))))
LOAR_fert_data <- LTREB_data_forfert %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"31"="1", "32"="2", "33"="3", "34"="4", "35"="5", "36"="6", "37"="7", "38"="8","39"="9", "40"="10"))))
POAL_fert_data <- LTREB_data_forfert %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"3"="1", "4"="2", "8"="3", "9"="4", "10"="5", "11"="6", "15"="7", "16"="8","17"="9", "19"="10","151"="11","152"="12","153"="13","154"="14","155"="15","156"="16","157"="17","158"="18"))))
POSY_fert_data <- LTREB_data_forfert %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"1"="1", "2"="2", "5"="3", "6"="4", "7"="5", "12"="6", "13"="7", "14"="8","18"="9", "20"="10","141"="11","142"="12","143"="13","144"="14","145"="15","146"="16","147"="17","148"="18", "149"="19", "150"="20"))))

# Create model matrices for each species year*endo and plot random effects
AGPE_yearendo_Xs_fert <- AGPE_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

AGPE_plot_Xs_fert <- AGPE_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELRI_yearendo_Xs_fert <- ELRI_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)


ELRI_plot_Xs_fert <- ELRI_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

ELVI_yearendo_Xs_fert <- ELVI_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t)%>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

ELVI_plot_Xs_fert <- ELVI_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


FESU_yearendo_Xs_fert <- FESU_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, year_endo2008_0 = 0, year_endo2008_1 = 0, .before = TRUE) %>% 
  add_column(year_endo2011_0 = 0, .before = c("year_endo2011_1"))

FESU_plot_Xs_fert <- FESU_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

LOAR_yearendo_Xs_fert <- LOAR_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

LOAR_plot_Xs_fert <- LOAR_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)

POAL_yearendo_Xs_fert <- POAL_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE) %>% 
  add_column(year_endo2013_0 = 0, .before = c("year_endo2013_1")) %>% 
  add_column(year_endo2014_0 = 0, .before = c("year_endo2014_1")) %>% 
  add_column(year_endo2015_0 = 0, year_endo2015_1 = 0, year_endo_2016_0 = 0, .before = c("year_endo2016_1"))

POAL_plot_Xs_fert <- POAL_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


POSY_yearendo_Xs_fert <- POSY_fert_data %>% 
  select(id, endo_01, year_t) %>% 
  unite("year_endo", year_t:endo_01, remove = FALSE) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = year_endo, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id, - endo_01, -year_t) %>% 
  add_column(year_endo2007_0 = 0, year_endo2007_1 = 0, .before = TRUE)

POSY_plot_Xs_fert <- POSY_fert_data %>% 
  select(id, plot_index) %>% 
  mutate(row_id = 1:n()) %>% 
  mutate(one = 1) %>% spread(key = plot_index, value = one, fill = 0, sep = "") %>% 
  ungroup() %>% select(-id, -row_id)


# Create data lists to be used for the Stan model
AGPE_fert_data_list <- list(flw_t = AGPE_fert_data$FLW_COUNT_T,
                            yearendo_Xs = AGPE_yearendo_Xs_fert,
                            plot_Xs = AGPE_plot_Xs_fert,
                            logsize_t = AGPE_fert_data$logsize_t,
                            origin_01 = AGPE_fert_data$origin_01,
                            endo_01 = AGPE_fert_data$endo_01,
                            endo_index = AGPE_fert_data$endo_index,
                            year_t = AGPE_fert_data$year_t_index,
                            plot = AGPE_fert_data$plot_index,
                            N = nrow(AGPE_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(AGPE_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(AGPE_fert_data$endo_01)))
str(AGPE_fert_data_list)

ELRI_fert_data_list <- list(flw_t = ELRI_fert_data$FLW_COUNT_T,
                            yearendo_Xs = ELRI_yearendo_Xs_fert,
                            plot_Xs = ELRI_plot_Xs_fert,
                            logsize_t = ELRI_fert_data$logsize_t,
                            origin_01 = ELRI_fert_data$origin_01,
                            endo_01 = ELRI_fert_data$endo_01,
                            endo_index = ELRI_fert_data$endo_index,
                            year_t = ELRI_fert_data$year_t_index,
                            plot = ELRI_fert_data$plot_index,
                            N = nrow(ELRI_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(ELRI_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(ELRI_fert_data$endo_01)))
str(ELRI_fert_data_list)

ELVI_fert_data_list <- list(flw_t = ELVI_fert_data$FLW_COUNT_T,
                            yearendo_Xs = ELVI_yearendo_Xs_fert,
                            plot_Xs = ELVI_plot_Xs_fert,
                            logsize_t = ELVI_fert_data$logsize_t,
                            origin_01 = ELVI_fert_data$origin_01,
                            endo_01 = ELVI_fert_data$endo_01,
                            endo_index = ELVI_fert_data$endo_index,
                            year_t = ELVI_fert_data$year_t_index,
                            plot = ELVI_fert_data$plot_index,
                            N = nrow(ELVI_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(ELVI_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(ELVI_fert_data$endo_01)))
str(ELVI_fert_data_list)

FESU_fert_data_list <- list(flw_t = FESU_fert_data$FLW_COUNT_T,
                            yearendo_Xs = FESU_yearendo_Xs_fert,
                            plot_Xs = FESU_plot_Xs_fert,
                            logsize_t = FESU_fert_data$logsize_t,
                            origin_01 = FESU_fert_data$origin_01,
                            endo_01 = FESU_fert_data$endo_01,
                            endo_index = FESU_fert_data$endo_index,
                            year_t = FESU_fert_data$year_t_index,
                            plot = FESU_fert_data$plot_index,
                            N = nrow(FESU_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(FESU_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(FESU_fert_data$endo_01)))
str(FESU_fert_data_list)

LOAR_fert_data_list <- list(flw_t = LOAR_fert_data$FLW_COUNT_T,
                            yearendo_Xs = LOAR_yearendo_Xs_fert,
                            plot_Xs = LOAR_plot_Xs_fert,
                            logsize_t = LOAR_fert_data$logsize_t,
                            origin_01 = LOAR_fert_data$origin_01,
                            endo_01 = LOAR_fert_data$endo_01,
                            endo_index = LOAR_fert_data$endo_index,
                            year_t = LOAR_fert_data$year_t_index,
                            plot = LOAR_fert_data$plot_index,
                            N = nrow(LOAR_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(LOAR_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(LOAR_fert_data$endo_01)))
str(LOAR_fert_data_list)
POAL_fert_data_list <- list(flw_t = POAL_fert_data$FLW_COUNT_T,
                            yearendo_Xs = POAL_yearendo_Xs_fert,
                            plot_Xs = POAL_plot_Xs_fert,
                            logsize_t = POAL_fert_data$logsize_t,
                            origin_01 = POAL_fert_data$origin_01,
                            endo_01 = POAL_fert_data$endo_01,
                            endo_index = POAL_fert_data$endo_index,
                            year_t = POAL_fert_data$year_t_index,
                            plot = POAL_fert_data$plot_index,
                            N = nrow(POAL_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(POAL_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 18L,
                            nEndo =   length(unique(POAL_fert_data$endo_01)))
str(POAL_fert_data_list)

POSY_fert_data_list <- list(flw_t = POSY_fert_data$FLW_COUNT_T,
                            yearendo_Xs = POSY_yearendo_Xs_fert,
                            plot_Xs = POSY_plot_Xs_fert,
                            logsize_t = POSY_fert_data$logsize_t,
                            origin_01 = POSY_fert_data$origin_01,
                            endo_01 = POSY_fert_data$endo_01,
                            endo_index = POSY_fert_data$endo_index,
                            year_t = POSY_fert_data$year_t_index,
                            plot = POSY_fert_data$plot_index,
                            N = nrow(POSY_fert_data),
                            K = 5L,
                            lowerlimit = as.integer(min(POSY_fert_data$FLW_COUNT_T)),
                            nYear = 11L,
                            nPlot = 20L,
                            nEndo =   length(unique(POSY_fert_data$endo_01)))
str(POSY_fert_data_list)

##############################################################################
####### Preparing datalists for Spikelet/inflorescence Kernel ------------------------------
##############################################################################
LTREB_data_for_spike <- LTREB_full %>%
  filter(!is.na(FLW_T)) %>% 
  filter(FLW_T>0)

dim(LTREB_data_for_seedmeans)

# Creating individual species data lists to be passed to the model

# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_spike_data <- LTREB_data_for_spike %>% 
  filter(species == "AGPE") %>% 
  filter(!is.na(SPIKEPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"111"="1", "112"="2", "113"="3", "114"="4", "115"="5", "116"="6", "117"="7", "118"="8","119"="9", "120"="10"))))


ELRI_spike_data <- LTREB_data_for_spike %>% #Elymus species are seed/inf
  filter(species == "ELRI") %>% 
  filter(!is.na(SEEDPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"101"="1", "102"="2", "103"="3", "104"="4", "105"="5", "106"="6", "107"="7", "108"="8","109"="9", "110"="10"))))

ELVI_spike_data <- LTREB_data_for_spike %>% #Elymus species are seed/inf
  filter(species == "ELVI") %>%
  filter(!is.na(SEEDPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"91"="1", "92"="2", "93"="3", "94"="4", "95"="5", "96"="6", "97"="7", "98"="8","99"="9", "100"="10"))))


FESU_spike_data <- LTREB_data_for_spike %>% 
  filter(species == "FESU") %>% 
  filter(!is.na(SPIKEPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"121"="1", "122"="2", "123"="3", "124"="4", "125"="5", "126"="6", "127"="7", "128"="8","129"="9", "130"="10"))))

LOAR_spike_data <- LTREB_data_for_spike %>% 
  filter(species == "LOAR") %>% 
  filter(!is.na(SPIKEPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"31"="1", "32"="2", "33"="3", "34"="4", "35"="5", "36"="6", "37"="7", "38"="8","39"="9", "40"="10"))))

POAL_spike_data <- LTREB_data_for_spike %>% 
  filter(species == "POAL") %>% 
  filter(!is.na(SPIKEPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"3"="1", "4"="2", "8"="3", "9"="4", "10"="5", "11"="6", "15"="7", "16"="8","17"="9", "19"="10","151"="11","152"="12","153"="13","154"="14","155"="15","156"="16","157"="17","158"="18"))))

POSY_spike_data <- LTREB_data_for_spike %>% 
  filter(species == "POSY") %>% 
  filter(!is.na(SPIKEPERINF_T)) %>% 
  mutate(plot_index = as.integer(as.character(recode(as.factor(plot_fixed),"1"="1", "2"="2", "5"="3", "6"="4", "7"="5", "12"="6", "13"="7", "14"="8","18"="9", "20"="10","141"="11","142"="12","143"="13","144"="14","145"="15","146"="16","147"="17","148"="18", "149"="19", "150"="20"))))

# Create data lists to be used for the Stan model
AGPE_spike_data_list <- list(spike_t = AGPE_spike_data$SPIKEPERINF_T,
                             logsize_t = AGPE_spike_data$logsize_t,
                             origin_01 = AGPE_spike_data$origin_01,
                             endo_01 = AGPE_spike_data$endo_01,
                             endo_index = AGPE_spike_data$endo_index,
                             year_t = AGPE_spike_data$year_t_index,
                             plot = AGPE_spike_data$plot_index,
                             N = length(na.omit(AGPE_spike_data$SPIKEPERINF_T)),
                             K = 5L,
                             nYear = 11L,
                             nPlot = 10L,
                             nEndo =   length(unique(AGPE_spike_data$endo_01)))
str(AGPE_spike_data_list)

ELRI_spike_data_list <- list(spike_t = (ELRI_spike_data$SEEDPERINF_T),
                             logsize_t = ELRI_spike_data$logsize_t,
                             origin_01 = ELRI_spike_data$origin_01,
                            endo_01 = ELRI_spike_data$endo_01,
                            endo_index = ELRI_spike_data$endo_index,
                            year_t = ELRI_spike_data$year_t_index,
                            plot = ELRI_spike_data$plot_index,
                            N = length((ELRI_spike_data$SEEDPERINF_T)),
                            N = nrow(ELRI_spike_data),
                            K = 5L,
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(ELRI_spike_data$endo_01)))
str(ELRI_spike_data_list)

ELVI_spike_data_list <- list(spike_t = (ELVI_spike_data$SEEDPERINF_T),
                             logsize_t = ELRI_spike_data$logsize_t,
                             origin_01 = ELRI_spike_data$origin_01,
                            endo_01 = ELVI_spike_data$endo_01,
                            endo_index = ELVI_spike_data$endo_index,
                            year_t = ELVI_spike_data$year_t_index,
                            plot = ELVI_spike_data$plot_index,
                            N = length((ELVI_spike_data$SEEDPERINF_T)),
                            N = nrow(ELVI_spike_data),
                            K = 5L,
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(ELVI_spike_data$endo_01)))
str(ELVI_spike_data_list)

FESU_spike_data_list <- list(spike_t = (FESU_spike_data$SPIKEPERINF_T),
                             logsize_t = FESU_spike_data$logsize_t,
                             origin_01 = FESU_spike_data$origin_01,
                            endo_01 = FESU_spike_data$endo_01,
                            endo_index = FESU_spike_data$endo_index,
                            year_t = FESU_spike_data$year_t_index,
                            plot = FESU_spike_data$plot_index,
                            N = length((FESU_spike_data$SPIKEPERINF_T)),
                            K = 5L,
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(FESU_spike_data$endo_01)))
str(FESU_spike_data_list)

LOAR_spike_data_list <- list(spike_t = (LOAR_spike_data$SPIKEPERINF_T),
                             logsize_t = LOAR_spike_data$logsize_t,
                             origin_01 = LOAR_spike_data$origin_01,
                            endo_01 = LOAR_spike_data$endo_01,
                            endo_index = LOAR_spike_data$endo_index,
                            year_t = LOAR_spike_data$year_t_index,
                            plot = LOAR_spike_data$plot_index,
                            N = length((LOAR_spike_data$SPIKEPERINF_T)),
                            K = 5L,
                            nYear = 11L,
                            nPlot = 10L,
                            nEndo =   length(unique(LOAR_spike_data$endo_01)))
str(LOAR_spike_data_list)

POAL_spike_data_list <- list(spike_t = (POAL_spike_data$SPIKEPERINF_T),
                             logsize_t = POAL_spike_data$logsize_t,
                             origin_01 = POAL_spike_data$origin_01,
                            endo_01 = POAL_spike_data$endo_01,
                            endo_index = POAL_spike_data$endo_index,
                            year_t = POAL_spike_data$year_t_index,
                            plot = POAL_spike_data$plot_index,
                            N = length((POAL_spike_data$SPIKEPERINF_T)),
                            K = 5L,
                            nYear = 11L,
                            nPlot = 18L,
                            nEndo =   length(unique(POAL_spike_data$endo_01)))
str(POAL_spike_data_list)

POSY_spike_data_list <- list(spike_t = (POSY_spike_data$SPIKEPERINF_T),
                             logsize_t = POSY_spike_data$logsize_t,
                             origin_01 = POSY_spike_data$origin_01,
                            endo_01 = POSY_spike_data$endo_01,
                            endo_index = POSY_spike_data$endo_index,
                            year_t = POSY_spike_data$year_t_index,
                            plot = POSY_spike_data$plot_index,
                            N = length(POSY_spike_data$SPIKEPERINF_T),
                            K = 5L,
                            nYear = 11L,
                            nPlot = 20L,
                            nEndo =   length(unique(POSY_spike_data$endo_01)))
str(POSY_spike_data_list)






##############################################################################
####### Preparing datalists for Seed Means Kernel ------------------------------
##############################################################################
LTREB_data_for_seedmeans <- LTREB_full %>%
  filter(!is.na(FLW_T)) %>% 
  filter(FLW_T>0)

dim(LTREB_data_for_seedmeans)
# View(LTREB_data_for_seedmeans)

# Creating individual species data lists to be passed to the model

# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_seed_data <- LTREB_data_for_seedmeans %>% 
  filter(species == "AGPE") %>% 
  filter(!is.na(SEEDPERSPIKE_T)) 
ELRI_seed_data <- LTREB_data_for_seedmeans %>% #Elymus species are seed/inf
  filter(species == "ELRI") %>% 
  filter(!is.na(SEEDPERINF_T))
ELVI_seed_data <- LTREB_data_for_seedmeans %>% #Elymus species are seed/inf
  filter(species == "ELVI") %>%
  filter(!is.na(SEEDPERINF_T))
FESU_seed_data <- LTREB_data_for_seedmeans %>% 
  filter(species == "FESU") %>% 
  filter(!is.na(SEEDPERSPIKE_T))
LOAR_seed_data <- LTREB_data_for_seedmeans %>% 
  filter(species == "LOAR") %>% 
  filter(!is.na(SEEDPERSPIKE_T))
POAL_seed_data <- LTREB_data_for_seedmeans %>% 
  filter(species == "POAL") %>% 
  filter(!is.na(SEEDPERSPIKE_T))
POSY_seed_data <- LTREB_data_for_seedmeans %>% 
  filter(species == "POSY") %>% 
  filter(!is.na(SEEDPERSPIKE_T))



# Create data lists to be used for the Stan model
AGPE_seed_data_list <- list(seed = AGPE_seed_data$SEEDPERSPIKE_T,
                            endo_01 = AGPE_seed_data$endo_01,
                            endo_index = AGPE_seed_data$endo_index,
                            year = AGPE_seed_data$year_t_index,
                            plot = AGPE_seed_data$plot_fixed,
                            Nseed = length(na.omit(AGPE_seed_data$SEEDPERSPIKE_T)),
                            K = 2,
                            nYear = length(unique(AGPE_seed_data$year_t_index)),
                            nPlot = length(unique(AGPE_seed_data$plot_fixed)),
                            nEndo =   length(unique(AGPE_seed_data$endo_01)))
str(AGPE_seed_data_list)

ELRI_seed_data_list <- list(seed = na.omit(ELRI_seed_data$SEEDPERINF_T),
                            endo_01 = na.omit(ELRI_seed_data$endo_01),
                            endo_index = ELRI_seed_data$endo_index,
                            year = ELRI_seed_data$year_t_index,
                            plot = ELRI_seed_data$plot_fixed,
                            Nseed = length(na.omit(ELRI_seed_data$SEEDPERINF_T)),
                            N = nrow(ELRI_seed_data),
                            K = 2,
                            nYear = length(unique(ELRI_seed_data$year_t_index)),
                            nPlot = length(unique(ELRI_seed_data$plot_fixed)),
                            nEndo =   length(unique(ELRI_seed_data$endo_01)))
str(ELRI_seed_data_list)

ELVI_seed_data_list <- list(seed = na.omit(ELVI_seed_data$SEEDPERINF_T),
                            endo_01 = na.omit(ELVI_seed_data$endo_01),
                            endo_index = ELVI_seed_data$endo_index,
                            year = ELVI_seed_data$year_t_index,
                            plot = ELVI_seed_data$plot_fixed,
                            Nseed = length(na.omit(ELVI_seed_data$SEEDPERINF_T)),
                            N = nrow(ELVI_seed_data),
                            K = 2,
                            nYear = length(unique(ELVI_seed_data$year_t_index)),
                            nPlot = length(unique(ELVI_seed_data$plot_fixed)),
                            nEndo =   length(unique(ELVI_seed_data$endo_01)))
str(ELVI_seed_data_list)

FESU_seed_data_list <- list(seed = na.omit(FESU_seed_data$SEEDPERSPIKE_T),
                            endo_01 = na.omit(FESU_seed_data$endo_01),
                            endo_index = FESU_seed_data$endo_index,
                            year = FESU_seed_data$year_t_index,
                            plot = FESU_seed_data$plot_fixed,
                            Nseed = length(na.omit(FESU_seed_data$SEEDPERSPIKE_T)),
                            K = 2,
                            nYear = length(unique(FESU_seed_data$year_t_index)),
                            nPlot = length(unique(FESU_seed_data$plot_fixed)),
                            nEndo =   length(unique(FESU_seed_data$endo_01)))
str(FESU_seed_data_list)

LOAR_seed_data_list <- list(seed = na.omit(LOAR_seed_data$SEEDPERSPIKE_T),
                            endo_01 = na.omit(LOAR_seed_data$endo_01),
                            endo_index = LOAR_seed_data$endo_index,
                            year = LOAR_seed_data$year_t_index,
                            plot = LOAR_seed_data$plot_fixed,
                            Nseed = length(na.omit(LOAR_seed_data$SEEDPERSPIKE_T)),
                            K = 2,
                            nYear = length(unique(LOAR_seed_data$year_t_index)),
                            nPlot = length(unique(LOAR_seed_data$plot_fixed)),
                            nEndo =   length(unique(LOAR_seed_data$endo_01)))
str(LOAR_seed_data_list)

POAL_seed_data_list <- list(seed = na.omit(POAL_seed_data$SEEDPERSPIKE_T),
                            endo_01 = na.omit(POAL_seed_data$endo_01),
                            endo_index = POAL_seed_data$endo_index,
                            year = POAL_seed_data$year_t_index,
                            plot = POAL_seed_data$plot_fixed,
                            Nseed = length(na.omit(POAL_seed_data$SEEDPERSPIKE_T)),
                            K = 2,
                            nYear = length(unique(POAL_seed_data$year_t_index)),
                            nPlot = length(unique(POAL_seed_data$plot_fixed)),
                            nEndo =   length(unique(POAL_seed_data$endo_01)))
str(POAL_seed_data_list)

POSY_seed_data_list <- list(seed = na.omit(POSY_seed_data$SEEDPERSPIKE_T),
                            endo_01 = na.omit(POSY_seed_data$endo_01),
                            endo_index = POSY_seed_data$endo_index,
                            year = POSY_seed_data$year_t_index,
                            plot = POSY_seed_data$plot_fixed,
                            Nseed = length(na.omit(POSY_seed_data$SEEDPERSPIKE_T)),
                            K = 2,
                            nYear = length(unique(POSY_seed_data$year_t_index)),
                            nPlot = length(unique(POSY_seed_data$plot_fixed)),
                            nEndo =   length(unique(POSY_seed_data$endo_01)))
str(POSY_seed_data_list)



##############################################################################
####### Preparing datalists for Seed to seedling Kernel ------------------------------
##############################################################################


