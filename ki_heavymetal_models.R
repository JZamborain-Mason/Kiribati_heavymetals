#modelling heavy metal concentrations in fish given heavy metal concentrations in algae
#spatially structured


# Load some libraries
library(dplyr)
library(tidyverse)
library(brms)
library(rstan)


#make rstan work with all cores
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())

#functions
panel_hist = function(x, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h = hist(x, plot = FALSE)
  breaks = h$breaks; nB = length(breaks)
  y = h$counts; y = y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
#Standardize function for continuous variables
standardise = function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#variance inflation factors  function for mixed models
vif.mer = function (fit) {
  ## adapted from rms::vif
  v <-vcov(fit)
  nam = names(fixef(fit))
  ## exclude intercepts
  ns = sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v = v[-(1:ns), -(1:ns), drop = FALSE]
    nam = nam[-(1:ns)]
  }
  d = diag(v)^0.5
  v = diag(solve(v/(d %o% d)))
  names(v) = nam
  v
}

#correlation function
panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y, method = "pearson",use = "complete.obs"))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
}

#load required data
algae = read.csv("XRFalgaefull_match.csv")
fish = read.csv("finalfishmetalfamily.csv")

#get mean algae concentrations per site (for sites that had multiple samples)
algae_mean<-algae %>% group_by(Site) %>%summarise(HABITAT=HABITAT[1],
                                                  ISLAND=ISLAND[1],
                                                  Hg.ppb_algae=mean(Hg.ppb,na.rm=T),
                                                  Cr.ppm_algae=mean(Cr.ppm,na.rm=T),
                                                  Fe.ppm_algae=mean(Fe.ppm,na.rm=T),
                                                  Co.ppm_algae=mean(Co.ppm,na.rm=T),
                                                  Ni.ppm_algae=mean(Ni.ppm,na.rm=T),
                                                  Cu.ppm_algae=mean(Cu.ppm,na.rm=T),
                                                  Zn.ppm_algae=mean(Zn.ppm,na.rm=T),
                                                  As.ppm_algae=mean(As.ppm,na.rm=T),
                                                  Se.ppm_algae=mean(Se.ppm,na.rm=T),
                                                  Cd.ppm_algae=mean(Cd.ppm,na.rm=T),
                                                  Pb.ppm_algae=mean(Pb.ppm,na.rm=T))

#checking correlations in algae
pairs(~., data=algae_mean %>%select(Hg.ppb_algae:Pb.ppm_algae),lower.panel=panel.cor )

#potentially combining sites
fish$combined_site<-paste(fish$ISLAND_CODE,fish$NEAR_SITE,sep="_")
unique(algae$Site)
unique(fish$MATCH_ALGAE)

#merge sites
fish <-merge(fish, algae_mean %>%select(Site,Hg.ppb_algae:Pb.ppm_algae), by.x="MATCH_ALGAE", by.y="Site",all.x=T)

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


#modelling fish concentrations

#muscle tissue____________________________________________________________
unique(fish$TISSUE_TYPE)
fish_mt<-fish %>% filter(TISSUE_TYPE=="VISC")
length(unique(fish_mt$SAMPLE_NUMBER))

#there are more samples than unique sampel numbers (Haley commented that this was done randomly whene xtracting metal concentrations
#due to some validation excercise)
#until we hear back on the best way to use these without introducing bias (e.g., some samples corrected otgers not),
#i will select one sample
fish_mt_unique <- fish_mt %>% group_by(SAMPLE_NUMBER) %>% slice_sample(n=1)

#get unique samples per species, genus, family
fish_mt_unique %>% group_by(FISH_NAME) %>% count()
fish_mt_unique %>% group_by(FISH_GENUS) %>% count()
fish_mt_unique %>% group_by(FAMILY) %>% count()
summary(as.factor(fish_mt_unique$DIET))
#standardize  potential covariates
fish_mt_unique$S_STA_LEG<-standardise(fish_mt_unique$STA_LEG)
fish_mt_unique$S_STA_LEG

#LETS TRY HG_______________________
summary(fish_mt_unique$Hg)
hist(fish_mt_unique$Hg)
hg_fish<-brm(Hg~Hg.ppb_algae+S_STA_LEG+DIET+(1|ISLAND_CODE/FAMILY/FISH_GENUS),data=fish_mt_unique[!is.na(fish_mt_unique$Hg.ppb_algae),],family=hurdle_lognormal())
pp_check(hg_fish)
plot(hg_fish)
conditional_effects(hg_fish)
#example effects
fixef(hg_fish)
ranef(hg_fish)
