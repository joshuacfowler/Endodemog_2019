## Authors: Josh and Tom	## Grass endophyte population model
## Purpose: Create a script that imports Endodemog data, perform all raw data manipulation to set up flowering tiller and seed production info,	
## and create an .RData object that can be loaded for analysis	
## Last Update: 03/12/2018
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)

# This is the main updated data sheet that we will merge with the flowering and seed data
LTREB_endodemog <- 
  read.csv(file = "endo_demog_long.csv")


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
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
  mutate(plot_fixed = (case_when(plot != "R" ~ plot, 
                                 plot == "R" ~ origin)))                       

# dim(LTREB_data)


##############################################################################
####### Preparing datalists for Survival Kernel ------------------------------
##############################################################################

# NA's in survival come from mostly 2017 recruits.
LTREB_data1 <- LTREB_data %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))

dim(LTREB_data1)

# LTREB_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + species)^2 + origin_01 
# , data = LTREB_data1)
# Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + species)^2 + origin_01 
# , data =LTREB_surv_for_matrix)

# 
# # Create data list for Stan model
# # # LTREB_surv_data_list <- list(surv_t1 = LTREB_data1$surv_t1,
# #                              Xs = Xs,    
# #                              endo_index = LTREB_data1$endo_index,
# #                              species_index = LTREB_data1$species_index,
# #                              spp_endo_index = LTREB_data1$spp_endo_index,
# #                              spp_year_index = LTREB_data1$spp_year_index,
# #                              year_t = LTREB_data1$year_t_index, 
# #                              N = nrow(LTREB_data1), 
# #                              K = ncol(Xs), 
# #                              nyear = length(unique(LTREB_data1$year_t_index)), 
# #                              nEndo =   length(unique(LTREB_data1$endo_01)),
# #                              nSpp = length(unique(LTREB_data1$species_index)))
# # /*** /*
# 
# str(LTREB_surv_data_list)
# 
# LTREB_sample <- sample_n(LTREB_data1, 1000)
# sample_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                                 , data = LTREB_sample)
# Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                    , data =sample_surv_for_matrix)
# 
# 
# # Create sample data list for Stan model
# sample_surv_data_list <- list(surv_t1 = LTREB_sample$surv_t1,
#                              Xs = Xs,    
#                              endo_index = LTREB_sample$endo_index,
#                              species_index = LTREB_sample$species_index,
#                              spp_endo_index = LTREB_sample$spp_endo_index,
#                              spp_year_index = LTREB_sample$spp_year_index,
#                              year_t = LTREB_sample$year_t_index, 
#                              N = nrow(LTREB_sample), 
#                              K = ncol(Xs), 
#                              nyear = length(unique(LTREB_sample$year_t_index)), 
#                              nEndo =   length(unique(LTREB_sample$endo_01)),
#                              nSpp = length(unique(LTREB_sample$species_index)))
# 
# 
# str(sample_surv_data_list)

# Creating individual species data lists to be passed to the model

# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_surv_data <- LTREB_data1 %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
ELRI_surv_data <- LTREB_data1 %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
ELVI_surv_data <- LTREB_data1 %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
FESU_surv_data <- LTREB_data1 %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
LOAR_surv_data <- LTREB_data1 %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
POAL_surv_data <- LTREB_data1 %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
POSY_surv_data <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))


# Create model matrices for each species
AGPE_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = AGPE_surv_data)
AGPE_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =AGPE_surv_for_matrix)

ELRI_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELRI_surv_data)
ELRI_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELRI_surv_for_matrix)

ELVI_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELVI_surv_data)
ELVI_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELVI_surv_for_matrix)

FESU_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = FESU_surv_data)
FESU_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =FESU_surv_for_matrix)

LOAR_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = LOAR_surv_data)
LOAR_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =LOAR_surv_for_matrix)

POAL_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POAL_surv_data)
POAL_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POAL_surv_for_matrix)

POSY_surv_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POSY_surv_data)
POSY_surv_matrix <- model.matrix(surv_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POSY_surv_for_matrix)

# Create data lists to be used for the Stan model
AGPE_surv_data_list <- list(surv_t1 = AGPE_surv_data$surv_t1,
                            Xs = AGPE_surv_matrix,
                            logsize_t = AGPE_surv_data$logsize_t,
                            origin_01 = AGPE_surv_data$origin_01,
                            endo_01 = AGPE_surv_data$endo_01,
                            endo_index = AGPE_surv_data$endo_index,
                            year_t = AGPE_surv_data$year_t_index,
                            plot = AGPE_surv_data$plot_index,
                            N = nrow(AGPE_surv_data),
                            K = ncol(AGPE_surv_matrix),
                            nYear = length(unique(AGPE_surv_data$year_t_index)),
                            nPlot = length(unique(AGPE_surv_data$plot_index)),
                            nEndo =   length(unique(AGPE_surv_data$endo_01)))
str(AGPE_surv_data_list)

ELRI_surv_data_list <- list(surv_t1 = ELRI_surv_data$surv_t1,
                            Xs = ELRI_surv_matrix,
                            logsize_t = ELRI_surv_data$logsize_t,
                            origin_01 = ELRI_surv_data$origin_01,
                            endo_01 = ELRI_surv_data$endo_01,
                            endo_index = ELRI_surv_data$endo_index,
                            year_t = ELRI_surv_data$year_t_index,
                            plot = ELRI_surv_data$plot_index,
                            N = nrow(ELRI_surv_data),
                            K = ncol(ELRI_surv_matrix),
                            nYear = length(unique(ELRI_surv_data$year_t_index)),
                            nPlot = length(unique(ELRI_surv_data$plot_index)),
                            nEndo =   length(unique(ELRI_surv_data$endo_01)))
str(ELRI_surv_data_list)

ELVI_surv_data_list <- list(surv_t1 = ELVI_surv_data$surv_t1,
                            Xs = ELVI_surv_matrix,
                            logsize_t = ELVI_surv_data$logsize_t,
                            origin_01 = ELVI_surv_data$origin_01,
                            endo_01 = ELVI_surv_data$endo_01,
                            endo_index = ELVI_surv_data$endo_index,
                            year_t = ELVI_surv_data$year_t_index,
                            plot = ELVI_surv_data$plot_index,
                            N = nrow(ELVI_surv_data),
                            K = ncol(ELVI_surv_matrix),
                            nYear = length(unique(ELVI_surv_data$year_t_index)),
                            nPlot = length(unique(ELVI_surv_data$plot_index)),
                            nEndo =   length(unique(ELVI_surv_data$endo_01)))
str(ELVI_surv_data_list)

FESU_surv_data_list <- list(surv_t1 = FESU_surv_data$surv_t1,
                            Xs = FESU_surv_matrix,
                            logsize_t = FESU_surv_data$logsize_t,
                            origin_01 = FESU_surv_data$origin_01,
                            endo_01 = FESU_surv_data$endo_01,
                            endo_index = FESU_surv_data$endo_index,
                            year_t = FESU_surv_data$year_t_index,
                            plot = FESU_surv_data$plot_index,
                            N = nrow(FESU_surv_data),
                            K = ncol(FESU_surv_matrix),
                            nYear = length(unique(FESU_surv_data$year_t_index)),
                            nPlot = length(unique(FESU_surv_data$plot_index)),
                            nEndo =   length(unique(FESU_surv_data$endo_01)))
str(FESU_surv_data_list)

LOAR_surv_data_list <- list(surv_t1 = LOAR_surv_data$surv_t1,
                            Xs = LOAR_surv_matrix,
                            logsize_t = LOAR_surv_data$logsize_t,
                            origin_01 = LOAR_surv_data$origin_01,
                            endo_01 = LOAR_surv_data$endo_01,
                            endo_index = LOAR_surv_data$endo_index,
                            year_t = LOAR_surv_data$year_t_index,
                            plot = LOAR_surv_data$plot_index,
                            N = nrow(LOAR_surv_data),
                            K = ncol(LOAR_surv_matrix),
                            nYear = length(unique(LOAR_surv_data$year_t_index)),
                            nPlot = length(unique(LOAR_surv_data$plot_index)),
                            nEndo =   length(unique(LOAR_surv_data$endo_01)))
str(LOAR_surv_data_list)

POAL_surv_data_list <- list(surv_t1 = POAL_surv_data$surv_t1,
                            Xs = POAL_surv_matrix,
                            logsize_t = POAL_surv_data$logsize_t,
                            origin_01 = POAL_surv_data$origin_01,
                            endo_01 = POAL_surv_data$endo_01,
                            endo_index = POAL_surv_data$endo_index,
                            year_t = POAL_surv_data$year_t_index,
                            plot = POAL_surv_data$plot_index,
                            N = nrow(POAL_surv_data),
                            K = ncol(POAL_surv_matrix),
                            nYear = length(unique(POAL_surv_data$year_t_index)),
                            nPlot = length(unique(POAL_surv_data$plot_index)),
                            nEndo =   length(unique(POAL_surv_data$endo_01)))
str(POAL_surv_data_list)

POSY_surv_data_list <- list(surv_t1 = POSY_surv_data$surv_t1,
                            Xs = POSY_surv_matrix,
                            logsize_t = POSY_surv_data$logsize_t,
                            origin_01 = POSY_surv_data$origin_01,
                            endo_01 = POSY_surv_data$endo_01,
                            endo_index = POSY_surv_data$endo_index,
                            year_t = POSY_surv_data$year_t_index,
                            plot = POSY_surv_data$plot_index,
                            N = nrow(POSY_surv_data),
                            K = ncol(POSY_surv_matrix),
                            nYear = length(unique(POSY_surv_data$year_t_index)),
                            nPlot = length(unique(POSY_surv_data$plot_index)),
                            nEndo =   length(unique(POSY_surv_data$endo_01)))
str(POSY_surv_data_list)





##############################################################################
####### Preparing datalists for Growth Kernel ------------------------------
##############################################################################

LTREB_data2 <- LTREB_data %>%
  filter(!is.na(logsize_t1)) %>%
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01)) %>% 
  mutate(size_t1 = as.integer(size_t1))

dim(LTREB_data2)
# LTREB_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                                 , data = LTREB_data2)
# Xs <- model.matrix(size_t1~ (logsize_t + endo_01 + species)^3 + origin_01 
#                    , data =LTREB_for_grow_matrix)
# 
# 
# # Create data list for Stan model
# LTREB_grow_data_list <- list(size_t1 = LTREB_data2$size_t1,
#                              Xs = Xs,    
#                              spp_endo_index = LTREB_data2$spp_endo_index,
#                              year_t = LTREB_data2$year_t_index, 
#                              N = nrow(LTREB_data2), 
#                              K = ncol(Xs), 
#                              lowerlimit = min(LTREB_data2$size_t1),
#                              nyear = length(unique(LTREB_data2$year_t_index)), 
#                              nEndo =   length(unique(LTREB_data2$endo_01)),
#                              nSpp = length(unique(LTREB_data2$species_index)))
# 
# 
# str(LTREB_grow_data_list)
# 
# LTREB_sample <- sample_n(LTREB_data2, 1000)
# sample_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                                  , data = LTREB_sample)
# Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                    , data =sample_for_grow_matrix)
# 
# 
# # Create sample data list for Stan model
# sample_grow_data_list <- list(size_t1 = LTREB_sample$size_t1,
#                               Xs = Xs,    
#                               spp_endo_index = LTREB_sample$spp_endo_index,
#                               year_t = LTREB_sample$year_t_index, 
#                               N = nrow(LTREB_sample), 
#                               K = ncol(Xs), 
#                               lowerlimit = min(LTREB_sample$size_t1),
#                               nyear = length(unique(LTREB_sample$year_t_index)), 
#                               nEndo =   length(unique(LTREB_sample$endo_01)),
#                               nSpp = length(unique(LTREB_sample$species_index)))
# 
# 
# str(sample_grow_data_list)


#########################################################################################################
# Creating individual species data lists to be passed to the model------------------------------
#########################################################################################################
# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_grow_data <- LTREB_data2 %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
ELRI_grow_data <- LTREB_data2 %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
ELVI_grow_data <- LTREB_data2 %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
FESU_grow_data <- LTREB_data2 %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
LOAR_grow_data <- LTREB_data2 %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
POAL_grow_data <- LTREB_data2 %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
POSY_grow_data <- LTREB_data2 %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))


# Create model matrices for each species
AGPE_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = AGPE_grow_data)
AGPE_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =AGPE_for_grow_matrix)

ELRI_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELRI_grow_data)
ELRI_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELRI_for_grow_matrix)

ELVI_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELVI_grow_data)
ELVI_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELVI_for_grow_matrix)

FESU_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = FESU_grow_data)
FESU_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =FESU_for_grow_matrix)

LOAR_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = LOAR_grow_data)
LOAR_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =LOAR_for_grow_matrix)

POAL_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POAL_grow_data)
POAL_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POAL_for_grow_matrix)

POSY_for_grow_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POSY_grow_data)
POSY_grow_matrix <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POSY_for_grow_matrix)

# Create data lists to be used for the Stan model
AGPE_grow_data_list <- list(size_t1 = AGPE_grow_data$size_t1,
                            Xs = AGPE_grow_matrix,
                            logsize_t = AGPE_grow_data$logsize_t,
                            origin_01 = AGPE_grow_data$origin_01,
                            endo_01 = AGPE_grow_data$endo_01,
                            endo_index = AGPE_grow_data$endo_index,
                            year_t = AGPE_grow_data$year_t_index,
                            plot = AGPE_grow_data$plot_index,
                            N = nrow(AGPE_grow_data),
                            K = ncol(AGPE_grow_matrix),
                            lowerlimit = as.integer(min(AGPE_grow_data$size_t1)),
                            nYear = length(unique(AGPE_grow_data$year_t_index)),
                            nPlot = length(unique(AGPE_grow_data$plot_index)),
                            nEndo =   length(unique(AGPE_grow_data$endo_01)))
str(AGPE_grow_data_list)

ELRI_grow_data_list <- list(size_t1 = ELRI_grow_data$size_t1,
                            Xs = ELRI_grow_matrix,
                            logsize_t = ELRI_grow_data$logsize_t,
                            origin_01 = ELRI_grow_data$origin_01,
                            endo_01 = ELRI_grow_data$endo_01,
                            endo_index = ELRI_grow_data$endo_index,
                            year_t = ELRI_grow_data$year_t_index,
                            plot = ELRI_grow_data$plot_index,
                            N = nrow(ELRI_grow_data),
                            K = ncol(ELRI_grow_matrix),
                            lowerlimit = as.integer(min(ELRI_grow_data$size_t1)),
                            nYear = length(unique(ELRI_grow_data$year_t_index)),
                            nPlot = length(unique(ELRI_grow_data$plot_index)),
                            nEndo =   length(unique(ELRI_grow_data$endo_01)))
str(ELRI_grow_data_list)

ELVI_grow_data_list <- list(size_t1 = ELVI_grow_data$size_t1,
                            Xs = ELVI_grow_matrix,
                            logsize_t = ELVI_grow_data$logsize_t,
                            origin_01 = ELVI_grow_data$origin_01,
                            endo_01 = ELVI_grow_data$endo_01,
                            endo_index = ELVI_grow_data$endo_index,
                            year_t = ELVI_grow_data$year_t_index,
                            plot = ELVI_grow_data$plot_index,
                            N = nrow(ELVI_grow_data),
                            K = ncol(ELVI_grow_matrix),
                            lowerlimit = as.integer(min(ELVI_grow_data$size_t1)),
                            nYear = length(unique(ELVI_grow_data$year_t_index)),
                            nPlot = length(unique(ELVI_grow_data$plot_index)),
                            nEndo =   length(unique(ELVI_grow_data$endo_01)))
str(ELVI_grow_data_list)

FESU_grow_data_list <- list(size_t1 = FESU_grow_data$size_t1,
                            Xs = FESU_grow_matrix,
                            logsize_t = FESU_grow_data$logsize_t,
                            origin_01 = FESU_grow_data$origin_01,
                            endo_01 = FESU_grow_data$endo_01,
                            endo_index = FESU_grow_data$endo_index,
                            year_t = FESU_grow_data$year_t_index,
                            plot = FESU_grow_data$plot_index,
                            N = nrow(FESU_grow_data),
                            K = ncol(FESU_grow_matrix),
                            lowerlimit = as.integer(min(FESU_grow_data$size_t1)),
                            nYear = length(unique(FESU_grow_data$year_t_index)),
                            nPlot = length(unique(FESU_grow_data$plot_index)),
                            nEndo =   length(unique(FESU_grow_data$endo_01)))
str(FESU_grow_data_list)

LOAR_grow_data_list <- list(size_t1 = LOAR_grow_data$size_t1,
                            Xs = LOAR_grow_matrix,
                            logsize_t = LOAR_grow_data$logsize_t,
                            origin_01 = LOAR_grow_data$origin_01,
                            endo_01 = LOAR_grow_data$endo_01,
                            endo_index = LOAR_grow_data$endo_index,
                            year_t = LOAR_grow_data$year_t_index,
                            plot = LOAR_grow_data$plot_index,
                            N = nrow(LOAR_grow_data),
                            K = ncol(LOAR_grow_matrix),            
                            lowerlimit = as.integer(min(LOAR_grow_data$size_t1)),
                            nYear = length(unique(LOAR_grow_data$year_t_index)),
                            nPlot = length(unique(LOAR_grow_data$plot_index)),
                            nEndo =   length(unique(LOAR_grow_data$endo_01)))
str(LOAR_grow_data_list)

POAL_grow_data_list <- list(size_t1 = POAL_grow_data$size_t1,
                            Xs = POAL_grow_matrix,
                            logsize_t = POAL_grow_data$logsize_t,
                            origin_01 = POAL_grow_data$origin_01,
                            endo_01 = POAL_grow_data$endo_01,
                            endo_index = POAL_grow_data$endo_index,
                            year_t = POAL_grow_data$year_t_index,
                            plot = POAL_grow_data$plot_index,
                            N = nrow(POAL_grow_data),
                            K = ncol(POAL_grow_matrix),
                            lowerlimit = as.integer(min(POAL_grow_data$size_t1)),
                            nYear = length(unique(POAL_grow_data$year_t_index)),
                            nPlot = length(unique(POAL_grow_data$plot_index)),
                            nEndo =   length(unique(POAL_grow_data$endo_01)))
str(POAL_grow_data_list)

POSY_grow_data_list <- list(size_t1 = POSY_grow_data$size_t1,
                            Xs = POSY_grow_matrix,
                            logsize_t = POSY_grow_data$logsize_t,
                            origin_01 = POSY_grow_data$origin_01,
                            endo_01 = POSY_grow_data$endo_01,
                            endo_index = POSY_grow_data$endo_index,
                            year_t = POSY_grow_data$year_t_index,
                            plot = POSY_grow_data$plot_index,
                            N = nrow(POSY_grow_data),
                            K = ncol(POSY_grow_matrix),
                            lowerlimit = as.integer(min(POSY_grow_data$size_t1)),
                            nYear = length(unique(POSY_grow_data$year_t_index)),
                            nPlot = length(unique(POSY_grow_data$plot_index)),
                            nEndo =   length(unique(POSY_grow_data$endo_01)))
str(POSY_grow_data_list)



########################################################################################################################
###### Reading in raw excel files to pull out reproduction data.---------------------------
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
LOAR_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_to2016_complete.xlsx", sheet = "LOAR")
LOAR_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_to2016_complete.xlsx", sheet = "LOAR recruits")
LOAR_data_r <- LOAR_data_r%>% # assign a Tag name for recruits and endophyte status based on plot and id
  mutate(tag = paste(Plot, `Recruit ID`, sep = "-")) %>% 
  mutate(Endo = LOAR_data$Endo[match(LOAR_data_r$Plot, LOAR_data$PLOT)])
LOAR_data_seed2008 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_to2016_complete.xlsx", sheet = "LOAR seeds 2008", skip=1)
LOAR_data_seed2009 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_to2016_complete.xlsx", sheet = "LOAR seeds 2009")

# Read in data from FESU
FESU_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU")
FESU_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU recruits")

# Read in data from ELVI
ELVI_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI(IN) originals up to 2016.xlsx", sheet = "ELVI")
ELVI_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVI recruits")
ELVI_data_seed <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVISeeds")
ELVI_data_r <- ELVI_data_r %>% 
  mutate(tag = paste(Plot, RecruitNo, sep = "-")) 
  

# Read in data from ELRI
ELRI_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI")
ELRI_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI recruits")
ELRI_data_seed2009 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI infl 09")
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




# Pulling out the seed production estimates. These are not measured for all plants, and so will go into a separate dataframe------------------------------
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


 # View(pseed)

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


pflw <- POAL_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(flw2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
pflw$year<- ifelse(pflw$variable == "flw2007", 2007, ifelse(pflw$variable == "Flwtillers1", 2008, ifelse(pflw$variable  == "FlwTillers2", 2009, ifelse(pflw$variable  == "FlwTillers3", 2010, ifelse(pflw$variable  == "FlwTillers4", 2011, ifelse(pflw$variable  == "FlwTillers5", 2012, ifelse(pflw$variable  == "FlwTillers6", 2013,ifelse(pflw$variable == "FlwTillers7", 2014,ifelse(pflw$variable == "FlwTillers8", 2015,ifelse(pflw$variable  == "FlwTillers9", 2016, NA))))))))))

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
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010","seed2011", "seed2012", "seed2013", "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))
rseed$year<- ifelse(rseed$variable == "seed2007", 2007, ifelse(rseed$variable == "seed2008", 2008, ifelse(rseed$variable == "seed2009", 2009, ifelse(rseed$variable == "seed2010", 2010, ifelse(rseed$variable  == "seed2011", 2011, ifelse(rseed$variable  == "seed2012", 2012, ifelse(rseed$variable  == "seed2013", 2013, ifelse(rseed$variable  == "seed2014", 2014, ifelse(rseed$variable  == "seed2015", 2015, ifelse(rseed$variable  == "seed2016", 2016, NA))))))))))

# View(rseed)

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
# View(rspike)

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


# View(pold_seed)

pold_spike <- POAL_data_old%>% 
  mutate("Birth Year" = year(`Date`)) %>% 
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
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))
pold_spike$year<- ifelse(pold_spike$variable == "spike2007", 2007, ifelse(pold_spike$variable == "Spikelets_InflA1", 2008, ifelse(pold_spike$variable == "Spikelets_InflB1", 2008, ifelse(pold_spike$variable == "spikelets_InflA3", 2009, ifelse(pold_spike$variable  == "spikelets_InflB3", 2009, ifelse(pold_spike$variable  == "spikelets_InflA4", 2010, ifelse(pold_spike$variable  == "spikelets_InflB4", 2010, ifelse(pold_spike$variable  == "spikelets_InflC4", 2010, ifelse(pold_spike$variable  == "spikelets_InflA5", 2011, ifelse(pold_spike$variable  == "spikelets_inflB5", 2011, ifelse(pold_spike$variable  == "spikelets_inflC5", 2011,ifelse(pold_spike$variable == "spikelets_inflA6", 2012, ifelse(pold_spike$variable == "spikelets_inflA7", 2013, ifelse(pold_spike$variable == "spikelets_InflB7", 2013, ifelse(pold_spike$variable == "spikelets_inflC7", 2013, ifelse(pold_spike$variable  == "spikelets_inflA8", 2014, ifelse(pold_spike$variable == "spikelets_InflB8", 2014, ifelse(pold_spike$variable == "spikelets_inflC8", 2014, ifelse(pold_spike$variable == "spikelets_inflA9", 2015, ifelse(pold_spike$variable == "spikelets_InflB9", 2015, ifelse(pold_spike$variable ==   "spikelets_inflC9", 2015, ifelse(pold_spike$variable == "spikelets_inflA10", 2016, ifelse(pold_spike$variable == "spikelets_InflB10", 2016, ifelse(pold_spike$variable == "spikelets_inflC10", 2016, NA))))))))))))))))))))))))
# View(pold_spike)

poldflw <- POAL_data_old %>% 
  mutate("Birth Year" = year(`Date`)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  mutate(flw2007 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
poldflw$year<- ifelse(poldflw$variable == "flw2007", 2007, ifelse(poldflw$variable == "Flwtillers1", 2008, ifelse(poldflw$variable  == "FlwTillers3", 2009, ifelse(poldflw$variable  == "FlwTillers4", 2010, ifelse(poldflw$variable  == "FlwTillers5", 2011, ifelse(poldflw$variable  == "FlwTillers6", 2012, ifelse(poldflw$variable  == "FlwTillers7", 2013,ifelse(poldflw$variable == "FlwTillers8", 2014,ifelse(poldflw$variable == "FlwTillers9", 2015,ifelse(poldflw$variable  == "FlwTillers10", 2016, NA))))))))))
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
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", 
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))
rold_seed$year<- ifelse(rold_seed$variable == "seed2007", 2007, ifelse(rold_seed$variable == "seed2008", 2008, ifelse(rold_seed$variable == "seed2008", 2008, ifelse(rold_seed$variable == "seed2009", 2009,ifelse(rold_seed$variable == "seed2010", 2010, ifelse(rold_seed$variable  == "seed2011", 2011, ifelse(rold_seed$variable  == "seed2012", 2012, ifelse(rold_seed$variable  == "seed2013", 2013, ifelse(rold_seed$variable  == "seed2014", 2014, ifelse(rold_seed$variable  == "seed2015", 2015, ifelse(rold_seed$variable  == "seed2016", 2016, NA)))))))))))

# View(rold_seed)

rold_spike <- POAL_data_old_r %>% 
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", "spikelets_inflA5", "spikelets_infl13","spike1","spike2", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("spike1",variable)~"A",
                              grepl("spike2", variable)~"B"))
rold_spike$year<- ifelse(rold_spike$variable == "spike2007", 2007, ifelse(rold_spike$variable == "spike2008", 2008, ifelse(rold_spike$variable == "spike2009", 2009, ifelse(rold_spike$variable == "spike2010", 2010, ifelse(rold_spike$variable == "spike2011", 2011, ifelse(rold_spike$variable == "spikelets_inflA5", 2012, ifelse(rold_spike$variable  == "spikelets_infl13", 2013, ifelse(rold_spike$variable  == "spike1", 2014, ifelse(rold_spike$variable == "spike2", 2014, ifelse(rold_spike$variable  == "spikelets_infl15", 2015, ifelse(rold_spike$variable  == "spikelets_infl16", 2016, NA)))))))))))
# View(rold_spike)

roldflw <- POAL_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = NA, flw2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "FLWTiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
roldflw$year<- ifelse(roldflw$variable == "flw2007", 2007, ifelse(roldflw$variable == "flw2008", 2008, ifelse(roldflw$variable == "FLWTiller09", 2009, ifelse(roldflw$variable == "FLWtiller10", 2010, ifelse(roldflw$variable == "FLWtiller11", 2011, ifelse(roldflw$variable  == "FLWtiller12", 2012, ifelse(roldflw$variable  == "FLWtiller13", 2013, ifelse(roldflw$variable  == "FLWtiller14", 2014, ifelse(roldflw$variable  == "FLWtiller15", 2015, ifelse(roldflw$variable  == "FLWtiller16", 2016, NA))))))))))
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


# View(po_seed)

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
# View(po_spike)


po_flw <- POSY_data %>% 
  mutate("Birth Year" = year(as.character(`Planted Date`))) %>% 
  mutate(flw2007 = NA, flw2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "flw2016"), 
       value.name = "flw") %>% 
  filter(!is.na(plot))
po_flw$year<- ifelse(po_flw$variable == "flw2007", 2007, ifelse(po_flw$variable == "Flwtillers1", 2008, ifelse(po_flw$variable  == "FlwTillers2", 2009, ifelse(po_flw$variable  == "FlwTillers3", 2010, ifelse(po_flw$variable  == "FlwTillers4", 2011, ifelse(po_flw$variable  == "FlwTillers5", 2012, ifelse(po_flw$variable  == "FlwTillers6", 2013,ifelse(po_flw$variable == "FlwTillers7", 2014,ifelse(po_flw$variable == "FlwTillers8", 2015,ifelse(po_flw$variable  == "flw2016", 2016, NA))))))))))
# View(po_flw)

po_seedmerge_ss <- left_join(po_spike, po_seed, by = c("plot", "pos", "tag", "Endo", 
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


# View(po_rseed)

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
# View(po_rspike)



po_rflw <- POSY_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_rflw$year<- ifelse(po_rflw$variable == "flw2007", 2007, ifelse(po_rflw$variable == "flw2008", 2008, ifelse(po_rflw$variable == "flw2009", 2009, ifelse(po_rflw$variable == "FLWtiller10", 2010, ifelse(po_rflw$variable == "FLWtiller11", 2011, ifelse(po_rflw$variable  == "FLWtiller12", 2012, ifelse(po_rflw$variable  == "FLWtiller13", 2013, ifelse(po_rflw$variable  == "FLWtiller14", 2014, ifelse(po_rflw$variable  == "FLWtiller15", 2015, ifelse(po_rflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(po_rflw)


po_rseedmerge_ss <- left_join(po_rspike, po_rseed, by = c( "plot", "pos", "tag", "Endo", 
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
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
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
# View(po_oldseed)

po_oldspike <- POSY_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
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
# View(po_oldspike)

po_oldflw <- POSY_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  mutate(flw2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n",
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
po_oldflw$year<- ifelse(po_oldflw$variable == "flw2007", 2007, ifelse(po_oldflw$variable == "FlwTillers1", 2008, ifelse(po_oldflw$variable  == "FlwTillers3", 2009, ifelse(po_oldflw$variable  == "FlwTillers4", 2010, ifelse(po_oldflw$variable  == "FlwTillers5", 2011, ifelse(po_oldflw$variable  == "FlwTillers6", 2012, ifelse(po_oldflw$variable  == "FlwTillers7", 2013,ifelse(po_oldflw$variable == "FlwTillers8", 2014,ifelse(po_oldflw$variable == "FlwTillers9", 2015,ifelse(po_oldflw$variable  == "FlwTillers10", 2016, NA))))))))))
# View(po_oldflw)


po_oldseedmerge_ss <- left_join(po_oldseed, po_oldspike, by = c("plot", "pos", "tag", "Endo", 
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
  mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                              grepl("spike2",variable) ~ "B",
                              grepl("spike3", variable) ~ "C"))
po_roldseed$year<- ifelse(po_roldseed$variable == "seed2007", 2007, ifelse(po_roldseed$variable == "seed2008", 2008, ifelse(po_roldseed$variable == "seed2009", 2009, ifelse(po_roldseed$variable == "seed2010", 2010, ifelse(po_roldseed$variable == "seed2011", 2011, ifelse(po_roldseed$variable  == "seed2012", 2012, ifelse(po_roldseed$variable  == "seed2013", 2013, ifelse(po_roldseed$variable  == "seed2014", 2014, ifelse(po_roldseed$variable  == "seed2015", 2015, ifelse(po_roldseed$variable  == "seed2016", 2016, NA))))))))))
# View(po_roldseed)


po_roldspike <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 =NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", 
                       "spike2011", "spikelets_inflA12","spikelets_inflA13",
                       "spikelets_inflb13","spikelets_inflb13",
                       "spike1", "spike2", "spike3",
                       "spike1__1", "spike2__1", "spike3__1",
                       "spike1__2", "spike2__2", "spike3__2"),
       value.name = "spikelets") %>% 
    mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                                grepl("spike2",variable) ~ "B",
                                grepl("spike3", variable) ~ "C"))
  
po_roldspike$year<- ifelse(po_roldspike$variable == "spike2007", 2007, ifelse(po_roldspike$variable == "spike2008", 2008, ifelse(po_roldspike$variable == "spike2009", 2009, ifelse(po_roldspike$variable == "spike2010", 2010, ifelse(po_roldspike$variable == "spike2011", 2011, ifelse(po_roldspike$variable  == "spikelets_inflA12", 2012, ifelse(po_roldspike$variable  == "spikelets_inflA13", 2013, ifelse(po_roldspike$variable  == "spikelets_inflb13", 2013, ifelse(po_roldspike$variable  == "spikelets_inflc13", 2013, ifelse(po_roldspike$variable  == "spike1", 2014, ifelse(po_roldspike$variable  == "spike2", 2014, ifelse(po_roldspike$variable  == "spike3", 2014, ifelse(po_roldspike$variable  == "spike1__1", 2015, ifelse(po_roldspike$variable  == "spike2__1", 2015, ifelse(po_roldspike$variable  == "spike3__1", 2015, ifelse(po_roldspike$variable  == "spike1__2", 2016, ifelse(po_roldspike$variable  == "spike2__2", 2016,ifelse(po_roldspike$variable  == "spike3__2", 2016,NA))))))))))))))))))
# View(po_roldspike)

po_roldflw <- POSY_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWTiller12", 
                       "FLWTiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_roldflw$year<- ifelse(po_roldflw$variable == "flw2007", 2007, ifelse(po_roldflw$variable == "flw2008", 2008, ifelse(po_roldflw$variable == "flw2009", 2009, ifelse(po_roldflw$variable == "FLWtiller10", 2010, ifelse(po_roldflw$variable == "FLWtiller11", 2011, ifelse(po_roldflw$variable  == "FLWTiller12", 2012, ifelse(po_roldflw$variable  == "FLWTiller13", 2013, ifelse(po_roldflw$variable  == "FLWtiller14", 2014, ifelse(po_roldflw$variable  == "FLWtiller15", 2015, ifelse(po_roldflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(po_roldflw)

po_roldseedmerge_ss <- left_join(po_roldspike, po_roldseed, by = c("plot", "pos", "tag", "Endo", 
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
LOAR_seed_tiller <- merge(LOAR_data, LOAR_data_seed2008, by.x = c("TAG"), by.y = c("Tag"), all.x = TRUE) #there is a separate sheet with the main raw data for seeds.
LOAR_seed_tiller1 <- merge(LOAR_seed_tiller, LOAR_data_seed2009, by.x = c("TAG"), by.y = c("Tag"), all.x = TRUE) #there is a separate sheet with the main raw data for seeds.


lseed <- LOAR_seed_tiller1 %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG" ) %>% 
  mutate(seed2007 = NA,
         seed2008a = Filled.x + totunfilled, seed2008b = Filled__1.x + totunfilled__1, 
         seed2009a = Filled.y + TOTunfilled, seed2009b = Filled__1.y + TOTunfilled__1,
         seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008a", "seed2008b", "seed2009a", "seed2009b",
                       "seed2010", "seed2011", "seed2012", "seed2013","seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
    mutate(tillerid = case_when(grepl("a", variable) ~ "A", 
                                grepl("b",variable) ~ "B",
                                grepl("c", variable) ~ "C"))
lseed$year<- ifelse(lseed$variable == "seed2007", 2007, ifelse(lseed$variable == "seed2008a", 2008, ifelse(lseed$variable == "seed2008b", 2008, ifelse(lseed$variable  == "seed2009a", 2009, ifelse(lseed$variable  == "seed2009b", 2009, ifelse(lseed$variable  == "seed2010", 2010, ifelse(lseed$variable  == "seed2011", 2011, ifelse(lseed$variable  == "seed2012", 2012, ifelse(lseed$variable  == "seed2013", 2013,ifelse(lseed$variable == "seed2014", 2014,ifelse(lseed$variable == "seed2015", 2015,ifelse(lseed$variable  == "seed2016", 2016, NA))))))))))))
# View(lseed)


lspike <- LOAR_seed_tiller1 %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(spike2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("spike2007", 
                       "Spikelets.x", "Spikelets__1.x",
                       "Spikelets.y", "Spikelets__1.y",
                       "spikelets_inflA3", "spikelets_inflB3",
                       "spikelets_inflA5", "spikelets_inflB5",
                       "spikelets_inflA6", "spikelets_inflB6",
                       "spikelets_inflA7", "spikelets_inflB7",
                       "spikelets_inflA8", "spikelets_inflB8",
                       "spikelets_inflA9", "spikelets_inflB9",
                       "spikelets_inflA10", "spikelets_inflB10"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("lets.", variable) ~ "A",
                              grepl("lets__1.", variable) ~ "B"))
lspike$year<- ifelse(lspike$variable == "spike2007", 2007, ifelse(lspike$variable == "Spikelets.x", 2008, ifelse(lspike$variable == "Spikelets__1.x", 2008, ifelse(lspike$variable == "Spikelets.y", 2009, ifelse(lspike$variable == "Spikelets__1.y", 2009, ifelse(lspike$variable  == "spike2", 2009, ifelse(lspike$variable  == "spikelets_inflA3", 2010, ifelse(lspike$variable  == "spikelets_inflB3", 2010, ifelse(lspike$variable  == "spikelets_inflA5", 2011, ifelse(lspike$variable  == "spikelets_inflB5", 2011, ifelse(lspike$variable  == "spikelets_inflA6", 2012, ifelse(lspike$variable  == "spikelets_inflB6", 2012, ifelse(lspike$variable  == "spikelets_inflA7", 2013, ifelse(lspike$variable  == "spikelets_inflB7", 2013, ifelse(lspike$variable == "spikelets_inflA8", 2014, ifelse(lspike$variable == "spikelets_inflB8", 2014, ifelse(lspike$variable == "spikelets_inflA9", 2015, ifelse(lspike$variable == "spikelets_inflB9", 2015, ifelse(lspike$variable  == "spikelets_inflA10", 2016, ifelse(lspike$variable  == "spikelets_inflB10", 2016, NA))))))))))))))))))))
# View(lspike)

lflw <- LOAR_seed_tiller1 %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(flw2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("flw2007", "FlwTillers1", "FLWTiller2", "FLWTillers3", 
                       "FLWTillers5", "FLWTillers6", "FLWTillers7", 
                       "FLWTillers8", "FLWTillers9", "FLWTillers10"), 
       value.name = "flw") 
lflw$year<- ifelse(lflw$variable == "flw2007", 2007, ifelse(lflw$variable == "FlwTillers1", 2008, ifelse(lflw$variable  == "FLWTiller2", 2009, ifelse(lflw$variable  == "FLWTillers3", 2010, ifelse(lflw$variable  == "FLWTillers5", 2011, ifelse(lflw$variable  == "FLWTillers6", 2012, ifelse(lflw$variable  == "FLWTillers7", 2013,ifelse(lflw$variable == "FLWTillers8", 2014,ifelse(lflw$variable == "FLWTillers9", 2015,ifelse(lflw$variable  == "FLWTillers10", 2016, NA))))))))))
# View(lflw)

l_seedmerge_ss <- left_join(lseed, lspike, by = c("plot", "pos", "tag", "Endo", 
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
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed")   %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
l_rseed$year<- ifelse(l_rseed$variable == "seed2007", 2007, ifelse(l_rseed$variable == "seed2008", 2008, ifelse(l_rseed$variable == "seed2009", 2009,ifelse(l_rseed$variable == "seed2010", 2010, ifelse(l_rseed$variable == "seed2011", 2011, ifelse(l_rseed$variable  == "seed2012", 2012, ifelse(l_rseed$variable  == "seed2013", 2013, ifelse(l_rseed$variable  == "seed2014", 2014, ifelse(l_rseed$variable  == "seed2015", 2015, ifelse(l_rseed$variable  == "seed2016", 2016, NA))))))))))
# View(l_rseed)

l_rspike <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>%
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
# View(l_rspike)

l_rflw <- LOAR_data_r %>%
  rename("Birth Year" = "Date", "plot" = "Plot", "pos" = "Recruit ID") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FlwTillers10","Flw11", "Flw12", 
                       "Flw13", "Flw14", "Flw15",
                       "Flw16"),
       value.name = "flw") 
l_rflw$year<- ifelse(l_rflw$variable == "flw2007", 2007, ifelse(l_rflw$variable == "flw2008", 2008, ifelse(l_rflw$variable == "flw2009", 2009, ifelse(l_rflw$variable == "FlwTillers10", 2010, ifelse(l_rflw$variable == "Flw11", 2011, ifelse(l_rflw$variable  == "Flw12", 2012, ifelse(l_rflw$variable  == "Flw13", 2013, ifelse(l_rflw$variable  == "Flw14", 2014, ifelse(l_rflw$variable  == "Flw15", 2015, ifelse(l_rflw$variable  == "Flw16", 2016, NA))))))))))
# View(l_rflw)

l_rseedmerge_ss <- left_join(l_rseed, l_rspike, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Birth Year", "year", "tillerid"))
# View(l_rseedmerge_ss)

l_rseedmerge_ssf <- merge(l_rseedmerge_ss, l_rflw, by = c( "plot", "pos", "tag", "Endo", 
                                                    "Birth Year", "year"),all = TRUE)
# View(l_rseedmerge_ssf)

# Combining recruit and original plant data
l_seedmerge_ssf <- po_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
                                       "tillerid", "flw","spikelets", "seed")]
l_rseedmerge_ssf <- po_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth Year","year",
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
# View(fseed)

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
# View(fspike)


fflw <- FESU_data %>% 
  rename("Birth Year" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = NA, flw2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "flw2008", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
fflw$year<- ifelse(fflw$variable == "flw2007", 2007, ifelse(fflw$variable == "flw2008", 2008, ifelse(fflw$variable  == "FlwTillers2", 2009, ifelse(fflw$variable  == "FlwTillers3", 2010, ifelse(fflw$variable  == "FlwTillers4", 2011, ifelse(fflw$variable  == "FlwTillers5", 2012, ifelse(fflw$variable  == "FlwTillers6", 2013,ifelse(fflw$variable == "FlwTillers7", 2014,ifelse(fflw$variable == "FlwTillers8", 2015,ifelse(fflw$variable  == "FlwTillers9", 2016, NA))))))))))
# View(fflw)

f_seedmerge_ss <- left_join(fseed, fspike, by = c("plot", "pos", "tag", "Endo", 
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
# View(f_rseed)

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
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

f_rspike$year<- ifelse(f_rspike$variable == "spike2007", 2007, ifelse(f_rspike$variable == "spike2008", 2008, ifelse(f_rspike$variable == "spike2009", 2009, ifelse(f_rspike$variable == "spike2010", 2010, ifelse(f_rspike$variable == "spike2011", 2011, ifelse(f_rspike$variable  == "spike2012", 2012, ifelse(f_rspike$variable  == "spikelets_infla13", 2013, ifelse(f_rspike$variable  == "spikelets_inflb13", 2013, ifelse(f_rspike$variable  == "spikelets_infla14", 2014, ifelse(f_rspike$variable  == "spikelets_inflb14", 2014, ifelse(f_rspike$variable  == "spikelets_infla15", 2015, ifelse(f_rspike$variable  == "spikelets_inflb15", 2015, ifelse(f_rspike$variable  == "spikelets_infla16", 2016, ifelse(f_rspike$variable  == "spikelets_inflb16", 2016, ifelse(f_rspike$variable  == "spikelets_inflc16", 2016,NA)))))))))))))))
# View(f_rspike)

f_rflw <- FESU_data_r %>%
  rename("Birth Year" = "Date", "Endo" = "endo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLW10", "FLW11", "FLW12", 
                       "FLW13", "FLW14", "FLW15",
                       "FLW16"),
       value.name = "flw") 
f_rflw$year<- ifelse(f_rflw$variable == "flw2007", 2007, ifelse(f_rflw$variable == "flw2008", 2008,ifelse(f_rflw$variable == "flw2009", 2009,ifelse(f_rflw$variable == "FLW10", 2010, ifelse(f_rflw$variable == "FLW11", 2011, ifelse(f_rflw$variable  == "FLW12", 2012, ifelse(f_rflw$variable  == "FLW13", 2013, ifelse(f_rflw$variable  == "FLW14", 2014, ifelse(f_rflw$variable  == "FLW15", 2015, ifelse(f_rflw$variable  == "FLW16", 2016, NA))))))))))
# View(f_rflw)

f_rseedmerge_ss <- left_join(f_rseed, f_rspike, by = c("plot", "pos", "tag", "Endo", 
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
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))

elviseed$year<- ifelse(elviseed$variable == "seed2007", 2007, ifelse(elviseed$variable == "Seeds_Infl08", 2008, ifelse(elviseed$variable  == "INfl1Seeds09", 2009, ifelse(elviseed$variable  == "Infl2Seeds09", 2009, ifelse(elviseed$variable  == "Infl3Seeds09", 2009, ifelse(elviseed$variable  == "Infl4Seeds09", 2009, ifelse(elviseed$variable  == "Infl1Seeds10", 2010, ifelse(elviseed$variable  == "Infl2Seeds10", 2010,ifelse(elviseed$variable  == "Infl1Seeds11", 2011,ifelse(elviseed$variable  == "Infl2Seeds11", 2011, ifelse(elviseed$variable  ==  "Seeds1_Infl12", 2012, ifelse(elviseed$variable  ==  "Seeds2_Infl12", 2012, ifelse(elviseed$variable  ==  "Seeds1_Infl13", 2013, ifelse(elviseed$variable  ==  "Seeds2_Infl13", 2013, ifelse(elviseed$variable =="Seeds1_Infl14", 2014, ifelse(elviseed$variable =="Seeds2_Infl14", 2014, ifelse(elviseed$variable == "Seeds1_Infl15", 2015, ifelse(elviseed$variable == "Seeds2_Infl15", 2015, ifelse(elviseed$variable  =="Seeds1_Infl16", 2016, ifelse(elviseed$variable  =="Seeds2_Infl16", 2016, NA))))))))))))))))))))
# View(elviseed)

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
# View(elvispike)

elviflw <- ELVI_seed_tiller %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = NA) %>% 
 melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers08", "FlwTillers09", "FlwTillers10", 
                       "FlwTillers11", "FlwTillers12", "FlwTillers13", 
                       "FlwTillers14", "FlwTillers15", "FlwTillers16"), 
       value.name = "flw") 
elviflw$year<- ifelse(elviflw$variable == "flw2007", 2007,ifelse(elviflw$variable == "FlwTillers08", 2008, ifelse(elviflw$variable  == "FlwTillers09", 2009, ifelse(elviflw$variable  == "FlwTillers10", 2010, ifelse(elviflw$variable  == "FlwTillers11", 2011, ifelse(elviflw$variable  == "FlwTillers12", 2012, ifelse(elviflw$variable  == "FlwTillers13", 2013,ifelse(elviflw$variable == "FlwTillers14", 2014,ifelse(elviflw$variable == "FlwTillers15", 2015,ifelse(elviflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elviflw)

elvi_seedmerge_ss <- left_join(elviseed, elvispike, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Birth Year","year","tillerid"))
# View(elvi_seedmerge_ss)

elvi_seedmerge_ssf <- merge(elvi_seedmerge_ss, elviflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth Year","year"), all = TRUE)
# View(elvi_merge_ssf)



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
                              grepl("fl2",variable) ~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))

elvi_rseed$year<- ifelse(elvi_rseed$variable == "seed2007", 2007,ifelse(elvi_rseed$variable == "seed2008", 2008,ifelse(elvi_rseed$variable == "seed2009", 2009, ifelse(elvi_rseed$variable == "seed2010", 2010, ifelse(elvi_rseed$variable == "seed2011", 2011, ifelse(elvi_rseed$variable  == "seed2012", 2012, ifelse(elvi_rseed$variable  == "seed2013", 2013, ifelse(elvi_rseed$variable  == "seed2014", 2014, ifelse(elvi_rseed$variable  == "seeds/infl1/15", 2015, ifelse(elvi_rseed$variable  == "seeds/infl2/15", 2015, ifelse(elvi_rseed$variable == "seed2016", 2016, NA)))))))))))
# View(elvi_rseed)

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
# View(elvi_rspike)

elvi_rflw <- ELVI_data_r %>%
  rename("Birth Year" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = NA, flw2008 =NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("flw2007", "flw2008", "FlwTillers09", "FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15"),
       value.name = "flw") 
elvi_rflw$year<- ifelse(elvi_rflw$variable == "flw2007", 2007, ifelse(elvi_rflw$variable == "flw2008", 2008, ifelse(elvi_rflw$variable == "FlwTillers09", 2009, ifelse(elvi_rflw$variable == "FlwTillers10", 2010, ifelse(elvi_rflw$variable == "FlwTillers11", 2011, ifelse(elvi_rflw$variable  == "FlwTillers12", 2012, ifelse(elvi_rflw$variable  == "FlwTillers13", 2013, ifelse(elvi_rflw$variable  == "FlwTillers14", 2014, ifelse(elvi_rflw$variable  == "FlwTillers15", 2015, ifelse(elvi_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elvi_rflw)

elvi_rseedmerge_ss <- left_join(elvi_rseed, elvi_rspike, by = c( "plot", "pos", "tag", "Endo", 
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
                              grepl("fl1",variable) ~ "A",
                              grepl("fl2",variable)~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))

elriseed$year<- ifelse(elriseed$variable == "seed2007", 2007, ifelse(elriseed$variable == "seed2008", 2008, ifelse(elriseed$variable  == "Infl1Seeds2", 2009, ifelse(elriseed$variable  == "Infl2Seeds2", 2009, ifelse(elriseed$variable  == "Infl3Seeds2", 2009, ifelse(elriseed$variable  == "Infl4Seeds2", 2009, ifelse(elriseed$variable  == "Seeds_Infl10", 2010, ifelse(elriseed$variable  == "seeds_1_11", 2011, ifelse(elriseed$variable  == "seeds_2_11", 2011, ifelse(elriseed$variable  == "seeds_1_12", 2012, ifelse(elriseed$variable  == "seeds_2_12", 2012,ifelse(elriseed$variable  == "seeds_1_13", 2013,ifelse(elriseed$variable  == "seeds_2_13", 2013, ifelse(elriseed$variable == "seeds_1_14", 2014,ifelse(elriseed$variable == "seeds_2_14", 2014, ifelse(elriseed$variable == "seeds_1_15", 2015, ifelse(elriseed$variable == "seeds_2_15", 2015,ifelse(elriseed$variable  == "Seeds1_Infl16", 2016,ifelse(elriseed$variable  == "Seeds2_Infl16", 2016, NA)))))))))))))))))))
# View(elriseed)

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
# View(elrispike)

elriflw <- ELRI_seed_tiller %>% 
  rename("Birth Year" = "Planted Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers08", "FlwTillers09.x", "FLwTillers10", 
                       "FlwTiller11", "FlwTiller12", "FlwTiller13", 
                       "FlwTiller14", "FlwTiller15", "FlwTillers16"), 
       value.name = "flw") 
elriflw$year<- ifelse(elriflw$variable == "flw2007", 2007, ifelse(elriflw$variable == "FlwTillers08", 2008, ifelse(elriflw$variable  == "FlwTillers09.x", 2009, ifelse(elriflw$variable  == "FLwTillers10", 2010, ifelse(elriflw$variable  == "FlwTiller11", 2011, ifelse(elriflw$variable  == "FlwTiller12", 2012, ifelse(elriflw$variable  == "FlwTiller13", 2013,ifelse(elriflw$variable == "FlwTiller14", 2014,ifelse(elriflw$variable == "FlwTiller15", 2015,ifelse(elriflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elriflw)

elri_seedmerge_ss <- left_join(elriseed, elrispike, by = c( "plot", "pos", "tag", "Endo", 
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
# View(elri_rseed)

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
# View(elri_rspike)

elri_rflw <- ELRI_data_r %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "FLWtiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLW13", "FLW14", "FLW15", "FlwTillers16"),
       value.name = "flw") 
elri_rflw$year<- ifelse(elri_rflw$variable == "seed2007", 2007, ifelse(elri_rflw$variable == "seed2008", 2008, ifelse(elri_rflw$variable == "FLWtiller09", 2009, ifelse(elri_rflw$variable == "FLWtiller10", 2010, ifelse(elri_rflw$variable == "FLWtiller11", 2011, ifelse(elri_rflw$variable  == "FLWtiller12", 2012, ifelse(elri_rflw$variable  == "FLW13", 2013, ifelse(elri_rflw$variable  == "FLW14", 2014, ifelse(elri_rflw$variable  == "FLW15", 2015, ifelse(elri_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elri_rflw)

elri_rseedmerge_ss <- left_join(elri_rseed, elri_rspike, by = c( "plot", "pos", "tag", "Endo", 
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
ELRIrepro$seed <- ifelse(ELRIrepro$seed == ".", NA, ELRIrepro$seed)

# ELRIrepro <- ELRIrepro[!(is.na(ELRIrepro$flw)),]

# View(ELRIrepro)







# Combining repro measurements across years for the AGPE data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken

agpeseed <- AGPE_data %>%
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011",  "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = NA)
agpeseed$year<-  ifelse(agpeseed$variable == "seed2007", 2007, ifelse(agpeseed$variable == "seed2008", 2008, ifelse(agpeseed$variable  == "seed2009", 2009, ifelse(agpeseed$variable  == "seed2010", 2010, ifelse(agpeseed$variable  == "seed2011", 2011, ifelse(agpeseed$variable  == "seed2012", 2012, ifelse(agpeseed$variable  == "seed2013", 2013,ifelse(agpeseed$variable == "seed2014", 2014,ifelse(agpeseed$variable == "seed2015", 2015,ifelse(agpeseed$variable  == "seed2016", 2016, NA))))))))))
# View(agpeseed)

agpespike <- AGPE_data %>% 
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012", "spike2013", 
                       "spike2014", "spike2015", "spike2016"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = NA)
agpespike$year<- ifelse(agpespike$variable == "spike2007", 2007, ifelse(agpespike$variable == "spike2008", 2008, ifelse(agpespike$variable  == "spike2009", 2009, ifelse(agpespike$variable  == "spike2010", 2010, ifelse(agpespike$variable  == "spike2011", 2011, ifelse(agpespike$variable  == "spike2012", 2012, ifelse(agpespike$variable  == "spike2013", 2013, ifelse(agpespike$variable == "spike2014", 2014, ifelse(agpespike$variable == "spike2015", 2015, ifelse(agpespike$variable  == "spike2016", 2016, NA))))))))))
# View(agpespike)

agpeflw <- AGPE_data %>% 
  rename("Birth Year" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth Year" = year(as.character(`Birth Year`))) %>% 
  mutate(flw2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FLWTiller3", 
                       "FLWTiller4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
agpeflw$year<- ifelse(agpeflw$variable == "flw2007", 2007, ifelse(agpeflw$variable == "Flwtillers1", 2008, ifelse(agpeflw$variable  == "FlwTillers2", 2009, ifelse(agpeflw$variable  == "FLWTiller3", 2010, ifelse(agpeflw$variable  == "FLWTiller4", 2011, ifelse(agpeflw$variable  == "FlwTillers5", 2012, ifelse(agpeflw$variable  == "FlwTillers6", 2013,ifelse(agpeflw$variable == "FlwTillers7", 2014,ifelse(agpeflw$variable == "FlwTillers8", 2015,ifelse(agpeflw$variable  == "FlwTillers9", 2016, NA))))))))))
# View(agpeflw)

agpe_seedmerge_ss <- merge(agpeseed, agpespike, by = c( "plot", "pos", "tag", "Endo", 
                                                    "Birth Year","year","tillerid"))
# View(agpe_seedmerge_ss)

agpe_seedmerge_ssf <- merge(agpe_seedmerge_ss, agpeflw, by = c("plot", "pos", "tag", "Endo", 
                                                       "Birth Year","year"), all = TRUE)
# View(agpe_seedmerge_ssf)


# Combining repro measurements across years for the AGPE recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
agpe_rseed <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed")  %>% 
  mutate(tillerid = NA)
agpe_rseed$year<- ifelse(agpe_rseed$variable == "seed2007", 2007, ifelse(agpe_rseed$variable == "seed2008", 2008, ifelse(agpe_rseed$variable == "seed2009", 2009, ifelse(agpe_rseed$variable == "seed2010", 2010, ifelse(agpe_rseed$variable == "seed2011", 2011, ifelse(agpe_rseed$variable  == "seed2012", 2012, ifelse(agpe_rseed$variable  == "seed2013", 2013, ifelse(agpe_rseed$variable  == "seed2014", 2014, ifelse(agpe_rseed$variable  == "seed2015", 2015, ifelse(agpe_rseed$variable == "seed2016", 2016, NA))))))))))
# View(agpe_rseed)

agpe_rspike <- AGPE_data_r %>%
  rename("Birth Year" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2013 = NA, spike2014 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("spike2007", "spike2008", "spike2009", 
                       "spike2010","spike2011", "SpikeletsA12", 
                       "spike2013", "spike2014", "spikepertillerA15","spikepertillerB15",
                       "spikepertillerA16","spikepertillerB16"),
       value.name = "spikelets")  %>% 
  mutate(tillerid = NA)
agpe_rspike$year<-  ifelse(agpe_rspike$variable == "spike2007", 2007, ifelse(agpe_rspike$variable == "spike2008", 2008, ifelse(agpe_rspike$variable == "spike2009", 2009, ifelse(agpe_rspike$variable == "spike2010", 2010, ifelse(agpe_rspike$variable == "spike2011", 2011, ifelse(agpe_rspike$variable  == "SpikeletsA12", 2012, ifelse(agpe_rspike$variable  == "spike2013", 2013, ifelse(agpe_rspike$variable  == "spike2014", 2014, ifelse(agpe_rspike$variable  == "spikepertillerA15", 2015, ifelse(agpe_rspike$variable  == "spikepertillerB15", 2015,ifelse(agpe_rspike$variable == "spikepertillerA16", 2016, ifelse(agpe_rspike$variable == "spikepertillerB16", 2016, NA))))))))))))
# View(agpe_rspike)

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

agpe_rseedmerge_ss <- left_join(agpe_rseed, agpe_rspike, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year", "tillerid"))
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



# AGPErepro <- AGPErepro[!(is.na(AGPErepro$flw)),]
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

# View(LTREB_repro)
LTREB_reprotem <- LTREB_repro %>% 
  group_by(plot, pos, tag, Endo, `Birth Year`, year, species, flw) %>% 
  summarize(mean_seed_per_tiller = mean(seed, na.rm = TRUE), mean_spike_per_tiller = mean(spikelets, na.rm = TRUE))


# View(LTREB_reprotem)

table(is.na(LTREB_repro_tem$seed), LTREB_reprotem$year, LTREB_repro$species)
table(LTREB_reprotem$year, is.na(LTREB_reprotem$mean_seed_per_tiller), LTREB_reprotem$species)
table(LTREB_repro$year, is.na(LTREB_repro$mean_spike_per_tiller), LTREB_repro$species)





## getting a dataframe with time t and t_1
LTREB_repro_t1 <-LTREB_repro %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, flwtillerno_t1 = flw, spikeperinfl_t1 = spikelets, seedperinfl_t1 = seed) %>%  
  mutate(year_t = year_t1 - 1)

LTREB_repro_t <- LTREB_repro %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, flwtillerno_t = flw, spikeperinfl_t = spikelets, seedperinfl_t = seed) %>% 
  mutate(year_t1 = year_t + 1)

LTREB_repro1 <- LTREB_repro_t1 %>% 
  full_join(LTREB_repro_t, by = c("plot", "pos", "tag", "Endo", "year_t", "year_t1", "species"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) 

# View(LTREB_repro1)
LTREB_temp <- merge(x = LTREB_data, y = LTREB_repro, by.x = c("plot", "pos","id", "endo", "birth", "year_t1", "species"), by.y = c("plot", "pos","tag", "Endo", "Birth Year", "year_t", "species"))






# getting flw tiller info to merge with endodemoglong
# agpeflw, agpe_rflw, elriflw, elri_rflw, elviflw, elvi_rflw,fflw,f_rflw, lflw, l_rflw, po_flw, po_rflw, po_oldflw, po_roldflw, pflw, poldflw, rflw, roldflw
agpeflw <- agpeflw %>% 
  mutate(species = "AGPE") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
agpe_rflw <- agpe_rflw %>% 
  mutate(species = "AGPE") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
elriflw <- elriflw%>% 
  mutate(species = "ELRI") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
elri_rflw <- elri_rflw%>% 
  mutate(species = "ELRI") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
elviflw <- elviflw%>% 
  mutate(species = "ELVI") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
elvi_rflw <- elvi_rflw%>% 
  mutate(species = "ELVI") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
fflw <- fflw%>% 
  mutate(species = "FESU") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
f_rflw <- f_rflw%>% 
  mutate(species = "FESU") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
lflw <- lflw%>% 
  mutate(species = "LOAR") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
l_rflw <- l_rflw%>% 
  mutate(species = "LOAR") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
po_flw <- po_flw%>% 
  mutate(species = "POSY") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
po_rflw <- po_rflw%>% 
  mutate(species = "POSY") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
po_oldflw <- po_oldflw%>% 
  mutate(species = "POSY") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
po_roldflw <- po_roldflw%>% 
  mutate(species = "POSY") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
pflw <- pflw%>% 
  mutate(species = "POAL") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
poldflw <- poldflw%>% 
  mutate(species = "POAL") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
rflw <- rflw%>% 
  mutate(species = "POAL") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)
roldflw <- roldflw%>% 
  mutate(species = "POAL") %>% 
  select(plot, pos, tag, Endo, flw, year, species, variable)


flwtiller <- agpeflw %>%  
  merge(agpe_rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(elriflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>%  
  merge(elri_rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(elviflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(elvi_rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>%  
  merge(fflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(f_rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(lflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(l_rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>%  
  merge(po_flw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(po_rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(po_oldflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(po_roldflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(pflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(poldflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) %>% 
  merge(rflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE)  %>% 
  merge(roldflw, by = c("plot", "pos","tag", "Endo", "flw", "year", "species", "variable"),all = TRUE) 

flwtiller1 <- flwtiller %>% 
  filter(is.na(flw) | flw != "?") %>% 
  filter(is.na(flw) | flw != "check tag")


## getting a dataframe with time t and t_1
flw_tiller_t1 <-flwtiller1 %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)

flw_tiller_t <- flwtiller1 %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, flw_t = flw) %>% 
  mutate(year_t1 = year_t + 1)

flw_tillermerge <- flw_tiller_t1 %>% 
  full_join(flw_tiller_t, by = c("plot", "pos", "tag", "Endo", "year_t", "year_t1", "species"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) 


# merge this with LTREB for flower model
LTREB_endodemog <- 
  read.csv(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
LTREB_tem <- LTREB_endodemog %>% 
  select(-seed_t, -seed_t1, -contains("spike")) %>% 
  rename("tag" = "id", "Endo" = "endo") %>% 
  mutate(Endo = case_when(Endo == "0" | Endo == "minus" ~ 0,
                             Endo == "1"| Endo =="plus" ~ 1))


# View(LTREB_tem)
LTREB_f <- merge(x = LTREB_tem, y = flw_tillermerge, by = c("plot", "pos", "tag", "species", "year_t1", "year_t", "Endo"), all = TRUE)

LTREB_f$flw_t1<- ifelse(is.na(LTREB_f$flw_t1.x) & is.na(LTREB_f$flw_t1.y), NA, ifelse(!is.na(LTREB_f$flw_t1.x) & is.na(LTREB_f$flw_t1.y), LTREB_f$flw_t1.x, ifelse(is.na(LTREB_f$flw_t1.x) & !is.na(LTREB_f$flw_t1.y), LTREB_f$flw_t1.y, NA)))

LTREB_flw <- LTREB_f %>% 
  mutate(flw_stat_t = case_when(flw_t == 0 ~0,
                                flw_t == 1 ~1)) %>% 
  mutate(flw_stat_t1 = case_when(flw_t1 == 0 ~ 0,
                                 flw_t1 > 0 ~ 1)) %>% 
  filter(!is.na(flw_stat_t1)) 
       
# View(LTREB_flw)
table(LTREB_flw$species, LTREB_flw$year_t)
table(LTREB_flw$flw_stat_t, LTREB_flw$year_t)
table(LTREB_flw$flw_stat_t1, LTREB_flw$year_t1)
