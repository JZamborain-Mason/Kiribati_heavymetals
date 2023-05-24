#modelling heavy metal concentrations in fish given heavy metal concentrations in algae
#spatially structured


# Load some libraries
library(dplyr)
library(tidyverse)
library(brms)


#make rstan work with all cores
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())


#load required data
algae = read.csv("XRFalgaefull_useme.csv")
fish = read.csv("finalfishmetalfamily.csv")

#potentially combining sites
fish$combined_site<-paste(fish$ISLAND_CODE,fish$NEAR_SITE,sep="_")
unique(algae$Site)
unique(paste(fish$ISLAND_CODE,fish$NEAR_SITE,sep="_"))

#summary sample size by island
algae %>%group_by(ISLAND) %>%summarise(n=n_distinct(Cd.ppm))
algae %>%group_by(ISLAND, HABITAT) %>%summarise(n=n_distinct(Cd.ppm))

#model example of island variability
cd_alg<-brm(Cd.ppm~HABITAT+(1|ISLAND),data=algae[!is.na(algae$Cd.ppm),],family=hurdle_lognormal())
pp_check(cd_alg)
ranef_cd_algae<-as.data.frame(ranef(cd_alg))
ranef_cd_algae$island<-row.names(ranef_cd_algae)
ggplot(ranef_cd_algae,aes(x=ISLAND.Estimate.Intercept,y=island,xmin= ISLAND.Q2.5.Intercept,xmax=ISLAND.Q97.5.Intercept))+geom_errorbar()+geom_point()+geom_vline(xintercept = 0,lty=2)+
  xlab("offset from global intercept (cd_ppm)")
