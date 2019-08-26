# Preliminaries -----------------------------------------------------------
## Key components: No longer sampling using MCMC update -- exhaustive evaulation of all potential GAM models. Adding in an option to allow for either using a fixed baseline OR using a new baseline and future for each of the potential GAM models

# Libraries
library(mgcv)  
library(tidyverse)
library(mvtnorm)
library(Smisc)
library(rgeos)
library(akima)
library(sf)
library(viridis)
library(cowplot)
library(corrr)
library(maptools)
library(ROCR)
library(hydroGOF)
library(sdm)
library(ggrepel)
library(ggmap)
library(raster)
library(tools)
library(randomForest)
library(broom)
library(snakecase)
library(openair)
library(randomcoloR)

# Data prep - Alaways run -------------------------------------------------
## Core functions
# NegativeLL function -- the likelihood of NEVA votes given the base and future, which comes out of a potential gam model. In this formulaion, we propose that the NEVA acts on the biomass model component. Not the presence component. 
loglike_bio_func<- function(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL){
  
  if(FALSE){
    gam.mod0.p = dat.use$mod.fitted.p[[1]]
    gam.mod0.b = dat.use$mod.fitted.b[[1]]
    season = season.use
    cand.params.b = cand.params.b
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0)
    nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0)
    nevaE = NULL
    dir.brks<- c(0.33, 0.67)
    fix.params<- fix.params
  }
  
  ilink.p<- family(gam.mod0.p)$linkinv
  ilink.b<- family(gam.mod0.b)$linkinv
  
  # NEVA conversion components
  # New baseline and future for specific candidate parameters
  base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
  fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  
  # Get predictions for presence model component -- no sampling
  base.p<- predict.gam(gam.mod0.p, newdata = base.preds.use, type = "link")
  fut.p<- predict.gam(gam.mod0.p, newdata = fut.preds.use, type = "link")
  
  # Get predictions for biomass model component, which does involve sampling
  lpmat.base.b<- predict.gam(gam.mod0.b, newdata = base.preds.use, type = "lpmatrix")
  base.bio<- as.numeric(lpmat.base.b %*% t(cand.params.b))
  lpmat.fut.b<- predict.gam(gam.mod0.b, newdata = fut.preds.use, type = "lpmatrix")
  fut.bio<- as.numeric(lpmat.fut.b %*% t(cand.params.b))
  
  # Combined predictions
  base.c<- as.numeric(base.p * base.bio)
  fut.c<- as.numeric(fut.p * fut.bio)
  #base.c<- as.numeric(base.bio)
  #fut.c<- as.numeric(fut.bio)
  
  # Likelihood of the GAM biomass model...uses normal distribution
  # For proposed GAM, get new baseline and future statistics
  bmn<- mean(base.c, na.rm = T)
  bsd<- sd(base.c, na.rm = T)
  fmn<- mean(fut.c, na.rm = T)
  fsd<- sd(fut.c, na.rm = T)
  
  # Calculate directional effect and vulnerability (sensitivity and exposure) breaks
  Dbrks<- quantile(na.omit(base.c), prob = dir.brks)
  Vbrks<- quantile(na.omit(base.c), prob = seq(from = 0, to = 1, length.out = 6))
  md<- length(Vbrks)/2+1
  
  ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
  if(!is.null(nevaD)){
    dtmp<- pnorm(Dbrks, fmn, fsd) # Cumulative distribution function, p(x) <= dbrks 0.33 quant and p(x) <= dbrks 0.67 quant from N(fmn, fsd) dist
    pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
    
    like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
  }
  
  ## Vulnerability piece: sensitivity and exposure
  if(!is.null(nevaS)){
    sentmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    sentmp<- c(sentmp[1], diff(sentmp), 1-sentmp[length(sentmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    send<- rep(NA, md)
    send[1]<- sentmp[md]
    
    for(k in 1:(md-1)){
      send[k+1]<- sentmp[md+k]+sentmp[md-k]
    }
    
    like.sens<- try(dmultinom(nevaS, prob = send, log = TRUE), silent = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
  }
  
  if(!is.null(nevaE)){
    exptmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    exptmp<- c(exptmp[1], diff(exptmp), 1-exptmp[length(exptmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    expd<- rep(NA, md)
    expd[1]<- exptmp[md]
    
    for(k in 1:(md-1)){
      expd[k+1]<- exptmp[md+k]+exptmp[md-k]
    }
    
    like.exp<- try(dmultinom(nevaE, prob = expd, log = TRUE), silent = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
  }
  
  # Calculate and return relevant likelihood
  if(is.null(nevaD)){
    if(class(like.sens) == "try-error" | class(like.exp) == "try-error") {
      likelihood<- NA
    } else {
      likelihood<- like.sens + like.exp
    }
  } else if(is.null(nevaS) & is.null(nevaE)){
    if(class(like.dir) == "try-error"){
      likelihood<- NA
    } else {
      likelihood<- like.dir
    }
  } else if(!is.null(nevaS) & is.null(nevaE)){
    classes<- c(class(like.dir), class(like.sens))
    if(any(classes == "try-error")){
      likelihood<- NA
    } else {
      likelihood<- like.dir + like.sens
    }
  } else {
    classes<- c(class(like.dir), class(like.sens), class(like.exp))
    if(any(classes == "try-error")){
      likelihood<- NA
    } else {
      likelihood<- like.dir + like.sens + like.exp
    }
  }
  return(likelihood)
}

# Probability of the biomass model 
logprior_bio_func<- function(gam.mod0.b, cand.params.b){
  # Likelihood of the GAM baseline future (fitted curves) for presence component
  # Gather coefficiences and varcov matrix, needed to define the normal dist of each candidate parameter value
  coefs.mod.b<- coef(gam.mod0.b)
  vcov.mod.b<- vcov(gam.mod0.b)
  
  # Likelihood of these candidate parameter values, given coef (mean) and sigma (se)
  cand.params.logprior.b<- dmvnorm(cand.params.b, mean = coefs.mod.b, sigma = vcov.mod.b, log = TRUE)
  return(data.frame("Bio.Prior" =  cand.params.logprior.b))
}

# Posterior (loglike + logprior)
logposterior_bio_func<- function(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL){
  # Likelihood of NEVA dir and vuln votes, given GAM baseline future (fitted curves)
  if(class(cand.params.b) == "numeric"){
    cand.params.b<- matrix(cand.params.b, nrow = 1, ncol = length(cand.params.b), byrow = T, dimnames = list(NULL, names(coef(gam.mod0.b))))
  }
  
  # Likelihood of NEVA votes given the biomass model
  loglike.out<- loglike_bio_func(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL)
  
  # Probability of prior (gam biomass model)
  logprior.out<- logprior_bio_func(gam.mod0.b, cand.params.b)
  
  # Add em all together
  return(data.frame("Likelihood" = loglike.out, "Bio.Prior" = logprior.out$Bio.Prior, "Posterior" = loglike.out + logprior.out$Bio.Prior))
}

# Bayes Algorithm
sdm_neva_bayes_bio<- function(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL){
  
  # Calculate likelihood, prior and posterior given candidate
  post<- logposterior_bio_func(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL)
  likes<- cbind(post$Likelihood, post$Bio.Prior, post$Posterior)
  return(likes)
}

temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

#### Loading in the data and some prep
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") 

# Additional processing/filtering -- remove Atlantic sturgeon as an endangered species
dat<- dat[dat$COMNAME != "ATLANTIC STURGEON", ]

# Observation summaries
dat.surv.yrs<- dat %>%
  filter(., PRESENCE > 0) %>%
  group_by(COMNAME, SEASON) %>%
  dplyr::summarize(., 
                   "Survey Years" = n_distinct(EST_YEAR))

spp.remove<- dat.surv.yrs[dat.surv.yrs$'Survey Years' < 23,]

dat.obs.peryear<- dat %>%
  filter(., PRESENCE > 0) %>%
  group_by(COMNAME, EST_YEAR, SEASON) %>%
  dplyr::summarize(., 
                   "Occurrence" = n_distinct(STATION)) %>%
  mutate(., "Occurrence Rate" = Occurrence/max(Occurrence)) %>%
  arrange(., desc(`Occurrence Rate`))

dat.obs.agg<- dat.obs.peryear %>%
  group_by(COMNAME, SEASON) %>%
  mutate(., "Avg Occurrence Rate" = mean(`Occurrence Rate`, na.rm = T)) %>%
  arrange(., desc(`Avg Occurrence Rate`))

dat.obs.peryear<- dat.obs.peryear %>%
  arrange(., desc(`Occurrence Rate`))


# Training vs. testing
train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither"))

# Bottom trawl strata
bstrat<- st_read("~/GitHub/COCA/Data/BottomTrawlStrata/BTS_Strata.shp")

# Get names of strata
bstrat.names<- unique(bstrat$STRATA)

# Reduce dataset to use tows within the bottom trawl strata
dat<- dat[dat$STRATUM %in% bstrat.names,]

# Training data by season: fall and then spring, scaling variables before fitting models and defining "presence" response for presence model component
dat.train.f<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) 

dat.train.f$PRESENCE.BIOMASS<- ifelse(dat.train.f$BIOMASS > 0, 1, 0)

# Need to keep mean and sd from rescale to use when we predict or project to other time periods
base.depth.mean.f<- mean(abs(dat.train.f$DEPTH))
base.depth.sd.f<- sd(abs(dat.train.f$DEPTH))
base.shelf.mean.f<- mean(dat.train.f$SHELF_POS)
base.shelf.sd.f<- sd(dat.train.f$SHELF_POS)
base.temp.mean.f<- mean(dat.train.f$SEASONALMU.OISST, na.rm = T)
base.temp.sd.f<- sd(dat.train.f$SEASONALMU.OISST, na.rm = T)
fall.rescale.df<- data.frame(SEASON = "FALL", mean.t = base.temp.mean.f, sd.t = base.temp.sd.f, mean.depth = base.depth.mean.f, sd.depth = base.depth.sd.f, mean.shelf = base.shelf.mean.f, sd.shelf = base.shelf.sd.f)

# Now spring...
dat.train.s<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) 

dat.train.s$PRESENCE.BIOMASS<- ifelse(dat.train.s$BIOMASS > 0, 1, 0)

# Temps to rescale other variables
base.depth.mean.sp<- mean(abs(dat.train.s$DEPTH))
base.depth.sd.sp<- sd(abs(dat.train.s$DEPTH))
base.shelf.mean.sp<- mean(dat.train.s$SHELF_POS)
base.shelf.sd.sp<- sd(dat.train.s$SHELF_POS)
base.temp.mean.sp<- mean(dat.train.s$SEASONALMU.OISST, na.rm = T)
base.temp.sd.sp<- sd(dat.train.s$SEASONALMU.OISST, na.rm = T)
spring.rescale.df<- data.frame(SEASON = "SPRING", mean.t = base.temp.mean.sp, sd.t = base.temp.sd.sp, mean.depth = base.depth.mean.sp, sd.depth = base.depth.sd.sp, mean.shelf = base.shelf.mean.sp, sd.shelf = base.shelf.sd.sp)

## Testing dataframes and applying correct scale to match rescaled variables used in the model fitting process
dat.test.f<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
dat.test.f$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
dat.test.f$SHELF_POS.Scale<- mapply(temp.scale, dat.test.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
dat.test.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.f$SEASONALMU.OISST, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

dat.test.f$PRESENCE.BIOMASS<- ifelse(dat.test.f$BIOMASS > 0, 1, 0)

## Testing dataframes
dat.test.s<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
dat.test.s$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.s$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
dat.test.s$SHELF_POS.Scale<- mapply(temp.scale, dat.test.s$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
dat.test.s$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.s$SEASONALMU.OISST, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

dat.test.s$PRESENCE.BIOMASS<- ifelse(dat.test.s$BIOMASS > 0, 1, 0)

# Create nested dataframes, one for testing, one for training
# Training
dat.train<- dat.train.f %>%
  bind_rows(., dat.train.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TRAIN.DATA") %>%
  arrange(COMNAME)

# Testing
dat.test<- dat.test.f %>%
  bind_rows(., dat.test.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TEST.DATA") %>%
  arrange(COMNAME) 

## We will also need correctly "rescaled" versions for the projection time periods
# Need base.preds, fut.preds, nevaD and nevaV
spring.preds = "~/GitHub/COCA/Data/spring.rast.preds06122018.rds"
spring.preds<- readRDS(spring.preds)

fall.preds = "~/GitHub/COCA/Data/fall.rast.preds06122018.rds"
fall.preds<- readRDS(fall.preds)

base.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.)))

base.preds.sp<- base.preds.sp %>%
  left_join(., spring.rescale.df, by = "SEASON")
base.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
base.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, base.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
base.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.sp$Baseline, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

base.preds.sp<- base.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) 
base.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
base.preds.f$SHELF_POS.Scale<- mapply(temp.scale, base.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
base.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.f$Baseline, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
base.preds.f<- base.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

## Future
fut.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Spring.2055.mu, Spring.2055.pct05, Spring.2055.pct95, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
fut.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
fut.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
fut.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055.mu, spring.rescale.df$mean.t, spring.rescale.df$sd.t)
fut.preds.sp$SEASONAL05.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055.pct05, spring.rescale.df$mean.t, spring.rescale.df$sd.t)
fut.preds.sp$SEASONAL95.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055.pct95, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

fut.preds.sp<- fut.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Fall.2055.mu, Fall.2055.pct05, Fall.2055.pct95, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
fut.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
fut.preds.f$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
fut.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055.mu, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
fut.preds.f$SEASONAL05.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055.pct05, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
fut.preds.f$SEASONAL95.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055.pct95, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

fut.preds.f<- fut.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

## Also going to want these as prediction dataframes to plot GAM smooths
# Prediction dataframe
pred.dat.f<- with(fut.preds$Data[[1]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("FALL", 1500)))
rescaled.dat.f<- with(fut.preds$Data[[1]],
                      data.frame("SST" = seq(min(Fall.2055.mu, na.rm = T), max(Fall.2055.mu, na.rm = T), length = 500),
                                 "Depth" = seq(min(abs(DEPTH), na.rm = T), max(abs(DEPTH), na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(SHELF_POS, na.rm = T), max(SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("FALL", 500)))

pred.dat.s<- with(fut.preds$Data[[2]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("SPRING", 1500)))
rescaled.dat.s<- with(fut.preds$Data[[2]],
                      data.frame("SST" = seq(min(Spring.2055.mu, na.rm = T), max(Spring.2055.mu, na.rm = T), length = 500),
                                 "Depth" = seq(min(abs(DEPTH), na.rm = T), max(abs(DEPTH), na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(SHELF_POS, na.rm = T), max(SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("SPRING", 500)))

pred.dat<- bind_rows(pred.dat.f, pred.dat.s)
rescaled.dat<- bind_rows(rescaled.dat.f, rescaled.dat.s)

# Bring in vulnerability assessment datasets
dir.dat<- read.csv("~/GitHub/COCA/Data/JHareDirectionalEffect.csv") %>%
  mutate(., "COMNAME" = toupper(Species)) %>%
  dplyr::select(., -Species) %>%
  group_by(COMNAME) %>%
  nest(., .key = "Dir.Data")
vuln.dat<- read.csv("~/GitHub/COCA/Data/JHareQualitativeDataResults.csv") %>%
  group_by(., Species, Attribute.Category) %>%
  summarise_at(., .vars = c("Low", "Moderate", "High", "Very.High"), mean) 
exp.dat<- vuln.dat %>%
  filter(., Attribute.Category == "Exposure.Factor") %>%
  mutate(., "COMNAME" = toupper(Species)) %>%
  data.frame(.) %>%
  dplyr::select(., -Species, -Attribute.Category) %>%
  group_by(COMNAME) %>%
  nest(., .key = "Exp.Data")
sens.dat<- vuln.dat %>%
  filter(., Attribute.Category == "Sensitivity.Attribute") %>%
  mutate(., "COMNAME" = toupper(Species)) %>%
  data.frame(.) %>%
  dplyr::select(., -Species, -Attribute.Category) %>% 
  group_by(COMNAME) %>%
  nest(., .key = "Sens.Data")
qual.dat<- dir.dat %>%
  left_join(., exp.dat, by = "COMNAME") %>%
  left_join(., sens.dat, by = "COMNAME")

# Bind to training data..
dat.full<- dat.train %>%
  left_join(., qual.dat, by = "COMNAME")

# Now to test data
dat.full<- dat.full %>%
  left_join(., dat.test, by = c("COMNAME", "SEASON"))


response = "Presence"
mod.fitted.p = dat.use$mod.fitted.p[[1]]
mod.fitted.b = dat.use$mod.fitted.b[[1]]
season = season.use
base.preds = base.preds
fut.preds = fut.preds
like.posts = like.posts.b
posts.samps = posts.df.b

ilink <- family(mod.fitted.b)$linkinv
# Get SDM predictions
# Get maximimum value and extract row from posts.samps

mod.ind.b<- like.posts.b$Iteration[which.max(like.posts.b$Value)]
best.fit.b<- posts.samps.b[mod.ind.b,]
best.fit.mat.b<- matrix(as.numeric(best.fit.b), nrow = 1, ncol = length(best.fit.b), byrow = T, dimnames = list(NULL, names(coef(dat.use$mod.fitted.b))))

# SDM
sdm.base<- base.preds$Data[[match(season, base.preds$SEASON)]]
sdm.base<- sdm.base %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
sdm.base<- predict.gam(mod.fitted.p, newdata = sdm.base, type = "response")


# Mean climate model SST
newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
newdat.mu<- newdat.mu %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
sdm.fut.mu<- predict.gam(mod.fitted.p, newdata = newdat.mu, type = "response")

# pct05 climate model
newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
newdat.05<- newdat.05 %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
sdm.fut.pct05<- predict.gam(mod.fitted.p, newdata = newdat.05, type = "response")

# pct95 climate model
newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
newdat.95<- newdat.95 %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
sdm.fut.pct95<- predict.gam(mod.fitted.p, newdata = newdat.95, type = "response")

# Getting combined biomass predictions
# Get maximimum value and extract row from posts.samps
mod.ind<- like.posts$Iteration[which.max(like.posts$Value)]
best.fit<- posts.samps[mod.ind,]
best.fit.mat<- matrix(as.numeric(best.fit), nrow = 1, ncol = length(best.fit), byrow = T, dimnames = list(NULL, names(coef(mod.fitted.b))))

# Make predictions with these values
lpmat.base<- predict.gam(mod.fitted.b, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "lpmatrix")
combo.map.base<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base * exp(ilink(as.numeric(lpmat.base %*% t(best.fit.mat)))))

# Mean climate model SST
newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
newdat.mu<- newdat.mu %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
lpmat.fut.mu<- predict.gam(mod.fitted.b, newdata = newdat.mu, type = "lpmatrix")
combo.map.fut.mu<- data.frame("x" = newdat.mu$x, "y" = newdat.mu$y, "pred" = sdm.fut.mu * exp(ilink(as.numeric(lpmat.fut.mu %*% t(best.fit.mat)))))

# pct05 climate model
newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
newdat.05<- newdat.05 %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
lpmat.fut.pct05<- predict.gam(mod.fitted.b, newdata = newdat.05, type = "lpmatrix")
combo.map.fut.pct05<- data.frame("x" = newdat.05$x, "y" = newdat.05$y, "pred.05" = sdm.fut.pct05 * exp(ilink(as.numeric(lpmat.fut.pct05 %*% t(best.fit.mat)))))

# pct95 climate model
newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
newdat.95<- newdat.95 %>%
  unnest() %>%
  dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
lpmat.fut.pct95<- predict.gam(mod.fitted.b, newdata = newdat.95, type = "lpmatrix")
combo.map.fut.pct95<- data.frame("x" = newdat.95$x, "y" = newdat.95$y, "pred.95" = sdm.fut.pct95 * exp(ilink(as.numeric(lpmat.fut.pct95 %*% t(best.fit.mat)))))

# Maps
# Spatial projections
proj.wgs84<- "+init=epsg:4326" #WGS84
proj.utm<- "+init=epsg:2960" #UTM 19

# NELME domaine
nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
st_crs(nelme)<- proj.wgs84
nelme.sp<- as(nelme, "Spatial")

#Bounds
xlim.use<- c(-77, -65)
ylim.use<- c(35, 45)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us <- raster::getData("GADM",country="USA",level=1)
us.states <- us[us$NAME_1 %in% states,]
us.states <- gSimplify(us.states, tol = 0.025, topologyPreserve = TRUE)
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]
ca.provinces <- gSimplify(ca.provinces, tol = 0.025, topologyPreserve = TRUE)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

# Ideally, want a bit more fine scale of a map for visual purposes. To do that, we will be doing some interpolation and then cropping the map to the NELME border. This can take a bit of time, so instead of doing it for EVERY dataset, lets do it once here so we can get a vector of the points to keep. Then just pass that into the plotting function.
coords.df<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y)
pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Combo
data.use<- combo.map.base
pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

# Limits
combo.range<- range(combo.map.base$pred, na.rm = T)
combo.fut.range<- range(combo.map.fut.mu$pred, na.rm = T)
combo.lwr.range<- range(combo.map.fut.pct05$pred.05, na.rm = T)
combo.upr.range<- range(combo.map.fut.pct95$pred.95, na.rm = T)
combo.lim<- c(round(min(combo.range,  combo.fut.range, combo.lwr.range, combo.upr.range), 2), round(max(combo.range, combo.fut.range, combo.lwr.range, combo.upr.range), 2))

# Continuous scale plot
plot.base.combo<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nBase", na.value = "white", limits = combo.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))

## Future
# Combo future
data.use<- combo.map.fut.mu
pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

plot.fut.combo<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nAverage future", na.value = "white", limits = combo.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))

# Lower
data.use<- combo.map.fut.pct05
pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

plot.fut.combo.lwr<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nColder future", na.value = "white", limits = combo.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))

# Upper
data.use<- combo.map.fut.pct95
pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

plot.fut.combo.upr<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nWarmer future", na.value = "white", limits = combo.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))

# Differences
# Getting limits
combo<- data.frame("x" = combo.map.fut.mu$x, "y" = combo.map.fut.mu$y, "pred" = combo.map.fut.mu$pred - combo.map.base$pred)
combo.range<- range(combo$pred, na.rm = T)
combo.lwr<- data.frame("x" = combo.map.fut.pct05$x, "y" = combo.map.fut.pct05$y, "pred" = combo.map.fut.pct05$pred - combo.map.base$pred)
combo.lwr.range<- range(combo.lwr$pred, na.rm = T)
combo.upr<- data.frame("x" = combo.map.fut.pct95$x, "y" = combo.map.fut.pct95$y, "pred" = combo.map.fut.pct95$pred - combo.map.base$pred)
combo.upr.range<- range(combo.upr$pred, na.rm = T)
diff.lim<- c(round(min(combo.range, combo.lwr.range, combo.upr.range), 2), round(max(combo.range, combo.lwr.range, combo.upr.range), 2))

# Combo
data.use<- combo
pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

plot.diff.combo<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_gradient2(name = "NEVA Biomass Diff\nAverage future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 

# Lower
data.use<- combo.lwr
pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

plot.diff.combo.lwr<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_gradient2(name = "NEVA Biomass Diff\nCold future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 

# Lower
data.use<- combo.upr
pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
row.names(coords.keep)<- NULL
pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
names(pred.df.use)<- c("X", "Y", "z")

plot.diff.combo.upr<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_gradient2(name = "NEVA Biomass Diff\nWarm future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 

out<- plot_grid(plot.base.combo, NULL, NULL, plot.fut.combo, plot.fut.combo.lwr, plot.fut.combo.upr, plot.diff.combo, plot.diff.combo.lwr, plot.diff.combo.upr, nrow = 3, align = "hv", scale = 1.1)

# Mapping -----------------------------------------------------------------
# Projections data set with mod.fitted objects
out.dir<- "~/GitHub/COCA/Results/NormalVoting_BiomassIncPresNoExposure_03152019/"
results<- read_rds(paste(out.dir, "SDMPredictions.rds", sep = ""))

# Set these for species stuff
plot.path<- "~/Desktop/SG_Lobster2019"
season.use<- "FALL"
species.use<- "AMERICAN LOBSTER"
xlim.use<- c(-72, -65)
ylim.use<- c(40.5, 45)

# Extract projections
proj.vals.base<- results %>%
  filter(., COMNAME == species.use, Proj.Class == "Baseline.combo.b" & SEASON == season.use) %>%
  unnest()
proj.diff<- results %>%
  filter(., COMNAME == species.use, Proj.Class == "Future_mean_diff.combo.b" & SEASON == season.use) %>%
  unnest()

# Rescaling for nicer maps...
# NELME domaine
nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
st_crs(nelme)<- proj.wgs84

# Interpolation grid
coords.df<- data.frame("x" = proj.vals.base$x, "y" = proj.vals.base$y)
pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Baseline
data.use<- proj.vals.base
pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$Projection))
pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp.base<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp.base))
row.names(coords.keep)<- NULL
pred.df.base<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp.base$z)))
names(pred.df.base)<- c("long", "lat", "z")

# Difference
data.use<- proj.diff
pred.df.diff<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$Projection))
pred.df.interp<- interp(pred.df.diff[,1], pred.df.diff[,2], pred.df.diff[,3], duplicate = "mean", extrap = TRUE,
                        xo=seq(-87.99457, -57.4307, length = 115),
                        yo=seq(22.27352, 48.11657, length = 133))
pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

# Clip to nelme
pred.df.temp.diff<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
coords.keep<- as.data.frame(st_coordinates(pred.df.temp.diff))
row.names(coords.keep)<- NULL
pred.df.diff<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp.diff$z)))
names(pred.df.diff)<- c("long", "lat", "z")

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us <- raster::getData("GADM",country="USA",level=1)
us.states <- us[us$NAME_1 %in% states,]
us.states <- gSimplify(us.states, tol = 0.025, topologyPreserve = TRUE)
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]
ca.provinces <- gSimplify(ca.provinces, tol = 0.025, topologyPreserve = TRUE)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)


plot.out.base<- ggplot() +
  # Here you'd make adjustments to the point colors...
  geom_tile(data = pred.df.base, aes(x = long, y = lat, fill = z)) +
  scale_fill_viridis(option = "viridis", name = "Biomass\nBase", na.value = "white") +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  #geom_tile(data = foot.pts, aes(x = x, y = y, color = layer), fill = NA, show.legend = FALSE, size= 0.75) +
  #scale_color_manual(name = "Fished cells", values = c(gmri.gray, NA)) +
  xlim(xlim.use) + 
  ylim(ylim.use) +
  coord_fixed(1.3) 

plot.out.diff<- ggplot() +
  # Here you'd make adjustments to the point colors...
  geom_tile(data = pred.df.diff, aes(x = long, y = lat, fill = z)) +
  scale_fill_gradient2(name = "Biomass Diff", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white") +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  #geom_tile(data = foot.pts, aes(x = x, y = y, color = layer), fill = NA, show.legend = FALSE, size= 0.75) +
  #scale_color_manual(name = "Fished cells", values = c(gmri.gray, NA)) +
  xlim(xlim.use) + 
  ylim(ylim.use) +
  coord_fixed(1.3) 

plot.out<- plot_grid(plot.out.base, plot.out.diff, nrow = 1, align = "hv")
ggsave(paste(plot.path, season.use, species.use, "maps.jpg", sep = ""), plot.out, width = 11, height = 8, units = "in")
