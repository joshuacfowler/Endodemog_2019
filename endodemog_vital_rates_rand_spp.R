library(tidyverse)

joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
#tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm634/Dropbox/EndodemogData/"

LTREB_full <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_full.csv"))
LTREB_data_forsurv <- LTREB_full %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))
