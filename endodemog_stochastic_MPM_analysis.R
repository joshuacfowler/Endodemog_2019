## purpose use the fitted vital rate models to compute the stochastic growth rate and
## decompose endophyte effects on mean versus variance
source("endodemog_stochastic_MPM_source.R")

# Read in data and stan objects -------------------------------------------
joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
#tompath <- "C:/Users/tm634/Dropbox/EndodemogData/"

# Read in data and stan objects -------------------------------------------
## big demography data frame
LTREB_full <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_full.csv"))
## seeds per spikelet data frame
LTREB_data_for_seedmeans <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_repro1.csv")) %>% 
  mutate(seed_per_spike = seed/spikelets) %>% 
  mutate(SEED_PER_SPIKE= case_when(species != "AGPE" ~ seed_per_spike,
                                   species == "AGPE" & tillerid_fixed == "multitillermean" ~ seed, # AGPE has some of its seeds data already recorded as seed/spike
                                   species == "AGPE" & tillerid_fixed != "multitillermean" ~ seed_per_spike)) %>% 
  filter(!is.na(SEED_PER_SPIKE))
## recruitment data frame


# stan fits for surv, grow, flow, fert, spikelets
surv_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/surv_fit_fixed.rds"))
grow_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/grow_fit_fixed.rds"))
flow_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/flow_fit_fixed.rds"))
fert_fit_endo_mean <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_fixed_endo_mean.rds"))
spike_fit_endo_mean <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/spike_fit_fixed_endo_mean.rds"))


surv_par <- extract(surv_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,tau_year))
grow_par <- extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,phi,tau_year))
flow_par <- extract(flow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,tau_year))

## try out the make_params function
## note that this is how species are indexing:
cbind(unique(LTREB_full$species),
      as.integer(as.numeric(as.factor(unique(LTREB_full$species)))))

## draw a posterior sample
i <- sample.int(10)

make_params(species,
            endo_mean,
            endo_var,
            original=0,
            rfx=F,
            year=NULL,
            surv_par,
            grow_par,
            flow_par,
            fert_par,
            spike_par,
            seed_par,
            recruit_par)
  






dim(surv_par$tau_year)
dim(flow_par$tau_year)

surv_par$beta0[100,]
surv_par$tau_year[3445,2,1,5]
surv_par$tau_year[3445,,,]

lapply(surv_par,mean)
