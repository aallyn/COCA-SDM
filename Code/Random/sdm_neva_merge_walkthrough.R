##### Walking through our SDM and NEVA merging process
#########

## Preliminary stuff
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

temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

## Function neccessities -- pieces needed for each of the functions
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")
spp.example<- c("ATLANTIC COD", "SUMMER FLOUNDER", "LONGFIN SQUID", "AMERICAN LOBSTER")
spp.example<- "ATLANTIC COD"

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") 

# Training vs. testing
train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither"))

# Bottom trawl strata
bstrat<- st_read("./Data/BottomTrawlStrata/BTS_Strata.shp")

# Get names of strata
bstrat.names<- unique(bstrat$STRATA)

# Reduce dataset
dat<- dat[dat$STRATUM %in% bstrat.names,]

# Training data
dat.train.f<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) %>%
  filter(., COMNAME %in% spp.example)
dat.train.f$PRESENCE.ABUNDANCE<- ifelse(dat.train.f$ABUNDANCE > 0, 1, 0)
dat.train.f$WT.ABUNDANCE<- ifelse(dat.train.f$PRESENCE.ABUNDANCE == 0, 1, dat.train.f$ABUNDANCE)

# Temps to rescale other time periods
base.depth.mean.f<- mean(abs(dat.train.f$DEPTH))
base.depth.sd.f<- sd(abs(dat.train.f$DEPTH))
base.shelf.mean.f<- mean(dat.train.f$SHELF_POS)
base.shelf.sd.f<- sd(dat.train.f$SHELF_POS)
base.temp.mean.f<- mean(dat.train.f$SEASONALMU.OISST, na.rm = T)
base.temp.sd.f<- sd(dat.train.f$SEASONALMU.OISST, na.rm = T)
fall.rescale.df<- data.frame(SEASON = "FALL", mean.t = base.temp.mean.f, sd.t = base.temp.sd.f, mean.depth = base.depth.mean.f, sd.depth = base.depth.sd.f, mean.shelf = base.shelf.mean.f, sd.shelf = base.shelf.sd.f)

dat.train.s<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) %>%
  filter(., COMNAME %in% spp.example)
dat.train.s$PRESENCE.ABUNDANCE<- ifelse(dat.train.s$ABUNDANCE > 0, 1, 0)
dat.train.s$WT.ABUNDANCE<- ifelse(dat.train.s$PRESENCE.ABUNDANCE == 0, 1, dat.train.s$ABUNDANCE)

# Temps to rescale other variables
base.depth.mean.sp<- mean(abs(dat.train.s$DEPTH))
base.depth.sd.sp<- sd(abs(dat.train.s$DEPTH))
base.shelf.mean.sp<- mean(dat.train.s$SHELF_POS)
base.shelf.sd.sp<- sd(dat.train.s$SHELF_POS)
base.temp.mean.sp<- mean(dat.train.s$SEASONALMU.OISST, na.rm = T)
base.temp.sd.sp<- sd(dat.train.s$SEASONALMU.OISST, na.rm = T)
spring.rescale.df<- data.frame(SEASON = "SPRING", mean.t = base.temp.mean.sp, sd.t = base.temp.sd.sp, mean.depth = base.depth.mean.sp, sd.depth = base.depth.sd.sp, mean.shelf = base.shelf.mean.sp, sd.shelf = base.shelf.sd.sp)

## Testing dataframes
dat.test.f<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
dat.test.f$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
dat.test.f$SHELF_POS.Scale<- mapply(temp.scale, dat.test.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
dat.test.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.f$SEASONALMU.OISST, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

dat.test.f<- dat.test.f %>%
  filter(., COMNAME %in% spp.example)
dat.test.f$PRESENCE.ABUNDANCE<- ifelse(dat.test.f$ABUNDANCE > 0, 1, 0)
dat.test.f$WT.ABUNDANCE<- ifelse(dat.test.f$PRESENCE.ABUNDANCE == 0, 1, dat.test.f$ABUNDANCE)

## Testing dataframes
dat.test.s<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
dat.test.s$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.s$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
dat.test.s$SHELF_POS.Scale<- mapply(temp.scale, dat.test.s$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
dat.test.s$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.s$SEASONALMU.OISST, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

dat.test.s<- dat.test.s %>%
  filter(., COMNAME %in% spp.example)
dat.test.s$PRESENCE.ABUNDANCE<- ifelse(dat.test.s$ABUNDANCE > 0, 1, 0)
dat.test.s$WT.ABUNDANCE<- ifelse(dat.test.s$PRESENCE.ABUNDANCE == 0, 1, dat.test.s$ABUNDANCE)

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

# Running the sdm_neva_bayes function. Easiest to map this function??
# Need base.preds, fut.preds, nevaD and nevaV
spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds03232018.rds"
spring.preds<- readRDS(spring.preds)

fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds03232018.rds"
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
  dplyr::select(., x, y, Spring.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
fut.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
fut.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
fut.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

fut.preds.sp<- fut.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Fall.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
fut.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
fut.preds.f$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
fut.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
fut.preds.f<- fut.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

## Fit the model
dat.train
dat.row<- 2
gam.mod0<- gam(PRESENCE.ABUNDANCE ~ s(DEPTH.Scale, k = 5, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, k = 5, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train$TRAIN.DATA[[dat.row]], family = binomial(link = logit), select = TRUE)
coef.mod0<- coef(gam.mod0)
vc.mod0<- vcov(gam.mod0)
season<- dat.train$SEASON[[dat.row]]

# Draws from posterior to get a mess of different candidate parameters
mod.sims<- 1000
cand.params.all<- rmvnorm(mod.sims, mean = coef.mod0, sigma = vc.mod0) 

## Some experimenting -- Same model just look at what happens when we change the voting...
base.update<- TRUE
dir.brks<- c(-1, 1)
vuln.brks<- c(-3, -2, -1, 1, 2, 3)

nevaD<- c(10, 2, 0)
nevaV<- c(1, 10, 1, 0)
dist<- "beta"
fix.params<- NULL

## Core functions
# Beta parameterization from mean and variance
est_beta_params_func<- function(mu, var) {
  alpha<- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta<- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# NegativeLL function -- the likelihood of NEVA votes given the base and future, which comes out of a potential gam model -- combined
loglike_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks, dist, fix.params){
  
  if(FALSE){
    gam.mod = gam.mod0
    season = season
    cand.params = matrix(cand.params.all[1,], nrow = 1, ncol = length(cand.params.all[1,]), byrow = T, dimnames = list(NULL, names(coef(gam.mod0)))) 
    cand.params = coef(gam.mod0)
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = nevaD
    nevaV = nevaV
    base.update = base.update
    dist = dist
    fix.params = fix.params
  }
  
  # NEVA conversion components
  # New baseline and future for specific candidate parameters
  # Get predictor matrix for each of the observations in base.preds and fut.preds
  base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
  fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  lpmat.base<- predict.gam(gam.mod0, newdata = base.preds.use, type = "lpmatrix")
  lpmat.fut<- predict.gam(gam.mod0, newdata = fut.preds.use, type = "lpmatrix")
  
  # Predictions, combining lpmats and candidate parameters
  # Multiply candidate parameters by predictor matrix
  if(base.update){
    base.prob<- lpmat.base %*% t(cand.params)
    fut.prob<- lpmat.fut %*% t(cand.params)
  } else {
    base.prob<- lpmat.base %*% coef(gam.mod0)
    fut.prob<- lpmat.fut %*% t(cand.params)
  }
  
  # Bin baseline distribution for directional effects and then bin vulnerability
  # Assuming normal distribution...
  if(dist == "norm"){
    # For proposed GAM, get new baseline and future statistics
    bmn<- mean(base.prob, na.rm = T)
    bsd<- sd(base.prob, na.rm = T)
    fmn<- mean(fut.prob, na.rm = T)
    fsd<- sd(fut.prob, na.rm = T)
    
    Dbrks<- dir.brks*bsd+bmn # Directional effect cuts (+/- 1SD), sets bins to integrate future pdf over
    Vbrks<- vuln.brks*bsd/2+bmn # Vulnerability cuts, sets bins to integrate future pdf over
    md<- length(Vbrks)/2+1 # Numeric value (2)
    
    ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
    if(!is.null(nevaD)){
      dtmp<- pnorm(Dbrks, fmn, fsd) # Cumulative distribution function, p(x) <= -1 SD and p(x) <= 1SD from N(fm2, fs2) dist
      pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
      
      like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    }
    
    ## Vulnerability piece
    if(!is.null(nevaV)){
      vtmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
      vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
      
      # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
      vd<- rep(NA, md)
      vd[1]<- vtmp[md]
      
      for(k in 1:(md-1)){
        vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      }
      
      like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    }
    
    if(is.null(nevaD)){
      likelihood<- like.vuln
    } else if(is.null(nevaV)){
      likelihood<- like.dir
    } else {
      likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
    }
  } 
  
  # Beta distribution
  if(dist == "beta"){
    ilink <- family(gam.mod0)$linkinv
    base.prob<- ilink(base.prob)
    fut.prob<- ilink(fut.prob)
    
    # For proposed GAM, get new baseline and future statistics
    bmn<- mean(base.prob, na.rm = T)
    bsd<- sd(base.prob, na.rm = T)
    fmn<- mean(fut.prob, na.rm = T)
    fsd<- sd(fut.prob, na.rm = T)
    
    # Parameterize beta distribution from mn and sd
    base.params<- est_beta_params_func(mu = bmn, var = bsd^2)
    fut.params<- est_beta_params_func(mu = fmn, var = fsd^2)
    
    Dbrks<- quantile(na.omit(base.prob), prob = c(0.37, 0.67))
    Vbrks<- quantile(na.omit(base.prob), prob = seq(from = 0, to = 1, length.out = 6))
    md<- length(Vbrks)/2+1
    
    ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
    if(!is.null(nevaD)){
      dtmp<- pbeta(Dbrks, fut.params$alpha, fut.params$beta)
      pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
      like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    }
    
    ## Vulnerability piece
    if(!is.null(nevaV)){
      vtmp<- pbeta(Vbrks, fut.params$alpha, fut.params$beta) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
      vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
      
      # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
      vd<- rep(NA, md)
      vd[1]<- vtmp[md]
      
      for(k in 1:(md-1)){
        vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      }
      
      like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    }
    
    if(is.null(nevaD)){
      likelihood<- like.vuln
    } else if(is.null(nevaV)){
      likelihood<- like.dir
    } else {
      likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
    }
  }
  
  if(dist == "norm.diff"){
    # For proposed GAM, get new baseline and future statistics
    diff<- base.prob - fut.prob
    diff.mn<- mean(diff, na.rm = T)
    diff.sd<- sd(diff, na.rm = T)
    
    # Quantiles
    Dbrks<- quantile(na.omit(diff), prob = c(0.25, 0.75))
    Vbrks<- quantile(na.omit(diff), prob = seq(from = 0, to = 1, length.out = 6))
    md<- length(Vbrks)/2+1
    
    ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
    if(!is.null(nevaD)){
      dtmp<- pnorm(Dbrks, diff.mn, diff.sd)
      pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
      like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    }
    
    ## Vulnerability piece
    if(!is.null(nevaV)){
      vtmp<- pnorm(Vbrks, diff.mn, diff.sd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
      vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
      
      # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
      vd<- rep(NA, md)
      vd[1]<- vtmp[md]
      
      for(k in 1:(md-1)){
        vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      }
      
      like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    }
    
    if(is.null(nevaD)){
      likelihood<- like.vuln
    } else if(is.null(nevaV)){
      likelihood<- like.dir
    } else {
      likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
    }
  }
  return(likelihood) # Loglikelihood
}


# Core internal functions are the loglike_func and logprior_func
loglike_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks){
  
  if(FALSE){
    gam.mod0 = gam.mod0
    season = season
    cand.params = coef.mod0 
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = nevaD
    nevaV = nevaV
    base.update = base.update
    dir.brks = dir.brks
    vuln.brks = vuln.brks
    testing = TRUE
  }
  
  # NEVA conversion components
  # New baseline and future for specific candidate parameters
  # Get predictor matrix for each of the observations in base.preds and fut.preds
  base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
  fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  lpmat.base<- predict.gam(gam.mod0, newdata = base.preds.use, type = "lpmatrix")
  lpmat.fut<- predict.gam(gam.mod0, newdata = fut.preds.use, type = "lpmatrix")
  
  # Predictions, combining lpmats and candidate parameters
  # Multiply candidate parameters by predictor matrix. Let's look at the differences here...
  cand.params<- matrix(cand.params, nrow = 1, ncol = length(coef(gam.mod0)), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
  if(testing){
    bases<- data.frame("bmn" = rep(NA, nrow(cand.params.all)), "bsd" = rep(NA, nrow(cand.params.all)))
    futs<- data.frame("fmn" = rep(NA, nrow(cand.params.all)), "fsd" = rep(NA, nrow(cand.params.all)))
    pd.out<- data.frame("Neg" = rep(NA, nrow(cand.params.all)), "Neut" = rep(NA, nrow(cand.params.all)), "Pos" = rep(NA, nrow(cand.params.all)))
    
    for(i in 1:nrow(cand.params.all)){
      cand.params.use<- matrix(cand.params.all[i,], nrow = 1, ncol = length(coef(gam.mod0)), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
      base.temp<- lpmat.base %*% t(cand.params.use)
      fut.temp<- lpmat.fut %*% t(cand.params.use)
      
      bases$bmn[i]<-  mean(base.temp, na.rm = T)
      bases$bsd[i]<-  sd(base.temp, na.rm = T)
      futs$fmn[i]<-  mean(fut.temp, na.rm = T)
      futs$fsd[i]<-  sd(fut.temp, na.rm = T)
      
      # Bin baseline distribution for directional effects and then bin vulnerability
      Dbrks<- dir.brks*bases$bsd[i]+bases$bmn[i] # Directional effect cuts (+/- 1SD), sets bins to integrate future pdf over
      Vbrks<- vuln.brks*bases$bsd[i]/2+bases$bmn[i] # Vulnerability cuts, sets bins to integrate future pdf over
      md<- length(Vbrks)/2+1 # Numeric value (2)
      
      dtmp<- pnorm(Dbrks, futs$fmn[i], futs$fsd[i]) # Cumulative distribution function, p(x) <= -1 SD and p(x) <= 1SD from N(fm2, fs2) dist
      pd.out[i,]<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
    }
    
    # Distribution of differences
    norm.diff<- data.frame("diff.mean" = bases$bmn - futs$fmn, "diff.sd" = bases$bsd + futs$fsd)
    hist(norm.diff$diff.mean)
    hist(rnorm(1000, mean = mean(norm.diff$diff.mean), sd = mean(norm.diff$diff.sd))) 
  }
  
  if(base.update){
    base.prob<- lpmat.base %*% t(cand.params)
    fut.prob<- lpmat.fut %*% t(cand.params)
  } else {
    base.prob<- lpmat.base %*% coef(gam.mod0)
    fut.prob<- lpmat.fut %*% t(cand.params)
  }
  
  # For proposed GAM, get new baseline and future statistics
  bmn<- mean(base.prob, na.rm = T)
  bsd<- sd(base.prob, na.rm = T)
  fmn<- mean(fut.prob, na.rm = T)
  fsd<- sd(fut.prob, na.rm = T)
  
  # Bin baseline distribution for directional effects and then bin vulnerability
  Dbrks<- dir.brks*bsd+bmn # Directional effect cuts (+/- 1SD), sets bins to integrate future pdf over
  Vbrks<- vuln.brks*bsd/2+bmn # Vulnerability cuts, sets bins to integrate future pdf over
  md<- length(Vbrks)/2+1 # Numeric value (2)
  
  ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
  if(!is.null(nevaD)){
    dtmp<- pnorm(Dbrks, fmn, fsd) # Cumulative distribution function, p(x) <= -1 SD and p(x) <= 1SD from N(fm2, fs2) dist
    pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
    
    like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
  }
  
  ## Vulnerability piece
  if(!is.null(nevaV)){
    vtmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    vd<- rep(NA, md)
    vd[1]<- vtmp[md]
    
    for(k in 1:(md-1)){
      vd[k+1]<- vtmp[md+k]+vtmp[md-k]
    }
    
    like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
  }
  
  if(is.null(nevaD)){
    likelihood<- like.vuln
  } else if(is.null(nevaV)){
    likelihood<- like.dir
  } else {
    likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
  }
  return(likelihood) # Loglikelihood
}

# Probability of our model 
logprior_func<- function(gam.mod0, cand.params, testing){
  if(FALSE){
    gam.mod0<- gam.mod0
    cand.params<- cand.params.all
    testing<- TRUE
  }
  # Likelihood of the GAM baseline future (fitted curves)
  # Gather coefficiences and varcov matrix, needed to define the normal dist of each candidate parameter value
  coefs.mod<- coef(gam.mod0) # Coefficients for every term/basis function in the model
  vcov.mod<- vcov(gam.mod0) 
  
  # Likelihood of each of the candidate parameter values, given coef (mean) and sigma (se)
  if(testing){
    mod.like<- rep(NA, nrow(cand.params.all))
    
    for(i in 1:nrow(cand.params.all)){
      cand.params.use<- matrix(cand.params.all[i,], nrow = 1, ncol = length(coef(gam.mod0)), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
      mod.like[i]<- dmvnorm(cand.params.use, mean = coefs.mod, sigma = vcov.mod, log = TRUE)
    }
    
    # When do these differences result in significant model differences??
    # Calculating the differences in deviance (-2*LL)
    temp<- data.frame(expand.grid(mod.like, mod.like))
    temp<- temp %>%
      rowwise() %>%
      mutate(., "ChiSq" = (-2*as.numeric(Var1))-(-2*as.numeric(Var2)),
             "PChiSq" = pchisq(ChiSq, df = 0, lower.tail = TRUE))
  }
  
  
  
  cand.params.logprior<- dmvnorm(cand.params, mean = coefs.mod, sigma = vcov.mod, log = TRUE)
  return(cand.params.logprior)
}
# Wrapper function is sdm_neva_bayes, this function calls the logposterior_function
# Bayes Algorithm
sdm_neva_bayes<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks){
  # Calculate likelihood, prior and posterior given candidate
  post<- logposterior_func(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks)
  likes<- cbind(post$Likelihood, post$Prior, post$Posterior)
  return(likes)
}

# Logposterior function -- this calls the loglike_func and the logprior_func
logposterior_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks){
  # Likelihood of NEVA dir and vuln votes, given GAM baseline future (fitted curves)
  if(class(cand.params) == "numeric"){
    cand.params<- matrix(cand.params, nrow = 1, ncol = length(cand.params), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
  }
  
  loglike.out<- loglike_func(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks)
  
  # Probability of prior
  logprior.out<- logprior_func(gam.mod0, cand.params)
  
  # Add em all together
  return(data.frame("Likelihood" = loglike.out, "Prior" = logprior.out, "Posterior" = loglike.out + logprior.out))
}



