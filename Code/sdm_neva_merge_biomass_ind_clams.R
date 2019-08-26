###### Merging SDM and NEVA: 07-06-2018
################
## Key components: No longer sampling using MCMC update -- exhaustive evaulation of all potential GAM models. Adding in an option to allow for either using a fixed baseline OR using a new baseline and future for each of the potential GAM models. 

## Preliminaries
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

#####
## Data prep
#####
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
dat.path<- "~/GitHub/COCA/Data/model.dat.clam.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  left_join(., fish.spp, by = "SVSPP") %>%
  filter(., season == "SUMMER") %>%
  arrange(., DATE)
names(dat)[names(dat) == "season"]<- "SEASON"
names(dat)<- toupper(names(dat))
dat$EST_YEAR<- format(as.Date(dat$DATE), "%Y")

# Add in sediment type...
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
neshelf.sediment<- readShapePoly("~/GitHub/COCA/Data/TNC_benthicsediment.shp") 
proj4string(neshelf.sediment)<- proj.wgs84

# Use over to get polygon data at each point
points.wgs84<- dat
coordinates(points.wgs84)<- ~DECDEG_BEGLON+DECDEG_BEGLAT
proj4string(points.wgs84)<- proj.wgs84
dat$SED.SIZE<- over(points.wgs84, neshelf.sediment)$SEDIMENT
dat$SED.TYPE<- over(points.wgs84, neshelf.sediment)$GRPSED

# Training vs. testing
train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither"))

# Training data: scaling variables before fitting models and defining "presence" response for presence model component
dat.train.s<- dat %>%
  filter(., TRAIN.TEST == "TRAIN") %>%
  mutate(., "YEAR" = factor(format(as.Date(DATE), "%Y"), levels = seq(from = min(format(as.Date(DATE), "%Y")), to = max(format(as.Date(DATE), "%Y")), by = 1)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) 

dat.train.s$PRESENCE.BIOMASS<- ifelse(dat.train.s$BIOMASS > 0, 1, 0)

# Need to keep mean and sd from rescale to use when we predict or project to other time periods
base.depth.mean.s<- mean(abs(dat.train.s$DEPTH))
base.depth.sd.s<- sd(abs(dat.train.s$DEPTH))
base.temp.mean.s<- mean(dat.train.s$SEASONALMU.OISST, na.rm = T)
base.temp.sd.s<- sd(dat.train.s$SEASONALMU.OISST, na.rm = T)
summ.rescale.df<- data.frame(SEASON = "SUMMER", mean.t = base.temp.mean.s, sd.t = base.temp.sd.s, mean.depth = base.depth.mean.s, sd.depth = base.depth.sd.s)


## Testing dataframes and applying correct scale to match rescaled variables used in the model fitting process
dat.test.s<- dat %>%
  filter(., TRAIN.TEST == "TEST") %>%
  mutate(., "YEAR" = factor(format(as.Date(DATE), "%Y"), levels = seq(from = min(format(as.Date(DATE), "%Y")), to = max(format(as.Date(DATE), "%Y")), by = 1))) %>%
  left_join(., summ.rescale.df, by = "SEASON")
dat.test.s$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.s$DEPTH), summ.rescale.df$mean.depth, summ.rescale.df$sd.depth)
dat.test.s$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.s$SEASONALMU.OISST, summ.rescale.df$mean.t, summ.rescale.df$sd.t)

dat.test.s$PRESENCE.BIOMASS<- ifelse(dat.test.s$BIOMASS > 0, 1, 0)

# Create nested dataframes, one for testing, one for training
# Training
dat.train<- dat.train.s %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TRAIN.DATA") %>%
  arrange(COMNAME)

# Testing
dat.test<- dat.test.s %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TEST.DATA") %>%
  arrange(COMNAME) 

## We will also need correctly "rescaled" versions for the projection time periods
# Need base.preds, fut.preds, nevaD and nevaV
summ.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/summer.rast.preds09122018.rds"
summ.preds<- readRDS(summ.preds)

base.preds.summ<- summ.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS, SED.SIZE, SED.TYPE) %>%
  mutate(., "SEASON" = rep("SUMMER", nrow(.)))

base.preds.summ<- base.preds.summ %>%
  left_join(., summ.rescale.df, by = "SEASON")
base.preds.summ$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.summ$DEPTH), summ.rescale.df$mean.depth, summ.rescale.df$sd.depth)
base.preds.summ$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.summ$Baseline, summ.rescale.df$mean.t, summ.rescale.df$sd.t)

base.preds.summ<- base.preds.summ %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.summ

## Future
fut.preds.summ<- summ.preds %>%
  dplyr::select(., x, y, Summer.2055.mu, Summer.2055.pct05, Summer.2055.pct95, DEPTH, SED.TYPE, SED.SIZE) %>%
  mutate(., "SEASON" = rep("SUMMER", nrow(.))) %>%
  left_join(., summ.rescale.df, by = "SEASON")
fut.preds.summ$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.summ$DEPTH), summ.rescale.df$mean.depth, summ.rescale.df$sd.depth)
fut.preds.summ$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.summ$Summer.2055.mu, summ.rescale.df$mean.t, summ.rescale.df$sd.t)
fut.preds.summ$SEASONAL05.OISST.Scale<- mapply(temp.scale, fut.preds.summ$Summer.2055.pct05, summ.rescale.df$mean.t, summ.rescale.df$sd.t)
fut.preds.summ$SEASONAL95.OISST.Scale<- mapply(temp.scale, fut.preds.summ$Summer.2055.pct95, summ.rescale.df$mean.t, summ.rescale.df$sd.t)

fut.preds.summ<- fut.preds.summ %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.summ 

## Also going to want these as prediction dataframes to plot GAM smooths
# Prediction dataframe
pred.dat.s<- with(fut.preds$Data[[1]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 500)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500)),
                             "SEASON" = rep("SUMMER", 1000)))
rescaled.dat.s<- with(fut.preds$Data[[1]],
                      data.frame("SST" = seq(min(Summer.2055.mu, na.rm = T), max(Summer.2055.mu, na.rm = T), length = 500),
                                 "Depth" = seq(min(abs(DEPTH), na.rm = T), max(abs(DEPTH), na.rm = T), length = 500),
                                 "Season" = rep("SUMMER", 1000)))

pred.dat<- bind_rows(pred.dat.s)
rescaled.dat<- bind_rows(rescaled.dat.s)

#####
## Model fit
#####
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

spp.model<- c("ATLANTIC SURFCLAM", "OCEAN QUAHOG", "NORTHERN QUAHOG")
dat.full<- dat.full %>%
  filter(., COMNAME %in% spp.model)

# Set up the loop
scenarios.type = "Dir.Sens"
mod.sims<- 1000
dir.brks<- c(0.33, 0.67)
out.dir.stem<- "~/Desktop/NormalVoting_BiomassIncPresNoExposure_01282019_SurfClam/"
fix.params<- NULL
set.seed(13)

if(dir.exists(out.dir.stem)){
  print("Directory Exists")
} else {
  dir.create(out.dir.stem)
  print("Directory Created")
}

gam_fit_full_func<- function(df, response){
  if(response == "Presence"){
    gam.mod0<- gam(PRESENCE.BIOMASS ~ SED.TYPE + s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = binomial(link = logit), select = TRUE)
    return(gam.mod0)
  } 
  
  if(response == "Biomass"){
    gam.mod0<- gam(BIOMASS.MOD ~ SED.TYPE + s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = gaussian, select = TRUE)
    return(gam.mod0)
  }
}

gam_fit_red_func<- function(df, response){
  if(response == "Presence"){
    gam.mod0<- gam(PRESENCE.BIOMASS ~ SED.TYPE + s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = binomial(link = logit), select = TRUE)
    return(gam.mod0)
  } 
  
  if(response == "Biomass"){
    gam.mod0<- gam(BIOMASS.MOD ~ SED.TYPE + s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = gaussian, select = TRUE)
    return(gam.mod0)
  }
}

resids_map_func<- function(test.data, response, predicted, type) {
  # Preliminary spatial stuff
  # Spatial projections
  proj.wgs84<- "+init=epsg:4326" #WGS84
  proj.utm<- "+init=epsg:2960" #UTM 19
  
  # NELME domaine
  nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
  st_crs(nelme)<- proj.wgs84
  
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
  
  if(response == "Presence" & type == "response"){
    resids<- predicted - test.data$PRESENCE.BIOMASS
    dat.resids<- data.frame("year" = test.data$EST_YEAR, "x" = test.data$DECDEG_BEGLON, "y" = test.data$DECDEG_BEGLAT, "resid" = resids)
    
    plots.out<- vector("list", length(unique(dat.resids$year)))
    
    zlim.use<- c(min(dat.resids$resid, na.rm = T), max(dat.resids$resid, na.rm = T))
    
    for(k in 1:length(unique(dat.resids$year))){
      data.use<- dat.resids[dat.resids$year == unique(dat.resids$year)[k],]
      pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$resid))
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
      
      # Plot
      plots.out[[k]]<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(unique(dat.resids$year)[k], " pred - obs", sep = ""), low = "blue", high = "red", mid = "white", midpoint = 0.0, na.value = "white", limits = c(-1,1)) +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "in"))
    }
    
    out<- plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], nrow = 2, ncol = 3, scale = 1)
  }
  
  if(response == "Biomass" & type == "response"){
    resids<- predicted - test.data$BIOMASS
    dat.resids<- data.frame("year" = test.data$EST_YEAR, "x" = test.data$DECDEG_BEGLON, "y" = test.data$DECDEG_BEGLAT, "resid" = resids)
    
    plots.out<- vector("list", length(unique(dat.resids$year)))
    
    zlim.use<- c(min(dat.resids$resid, na.rm = T), max(dat.resids$resid, na.rm = T))
    
    for(k in 1:length(unique(dat.resids$year))){
      data.use<- dat.resids[dat.resids$year == unique(dat.resids$year)[k],]
      pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$resid))
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
      
      # Plot
      plots.out[[k]]<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(unique(dat.resids$year)[k], " pred - obs", sep = ""), low = "blue", high = "red", mid = "white", midpoint = 0.0, na.value = "white") +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "in"))
    }
    
    out<- plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], nrow = 2, ncol = 3, scale = 1)
  }
  
  return(out)
}

predict_func<- function(mod.fitted.p, mod.fitted.b, response, percentile, test.data) {
  temp<- dplyr::select(test.data, one_of(c("SED.TYPE", "DEPTH.Scale", percentile)))
  test.data<- data.frame(na.omit(temp))
  out.p<- round(as.numeric(predict.gam(mod.fitted.p, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3)
  
  if(response == "Biomass"){
    out.b<- exp(round(as.numeric(predict.gam(mod.fitted.b, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3))
    out.c<- out.p * out.b
    return(out.c)
  } else {
    return(out.p)
  }
}

neva_predict_bio_func<- function(mod.fitted.p, mod.fitted.b, test.data, neva.best.fit.b = best.fit.mat.b){
  temp<- dplyr::select(test.data, one_of(c("DECDEG_BEGLON", "DECDEG_BEGLAT", "SED.TYPE", "DEPTH.Scale", "SEASONALMU.OISST.Scale")))
  temp.data<- data.frame(na.omit(temp))
  
  # Presence prediction
  neva.pred.p<- predict.gam(mod.fitted.p, newdata = temp.data, type = "response")
  
  # Biomass prediction
  lpmat<- predict.gam(mod.fitted.b, newdata = temp.data, type = "lpmatrix")
  neva.pred.b<- as.numeric(lpmat %*% t(neva.best.fit.b))
  ilink<- family(mod.fitted.b)$linkinv
  neva.pred.b<- exp(ilink(neva.pred.b)) 
  return(neva.pred.p*neva.pred.b)
}

pred_ranges_func<- function(predicted){
  pred.ranges<- data.frame("Min.Pred" = min(predicted, na.rm = T), "Max.Pred" = max(predicted, na.rm = T), "Mean.Pred" = mean(predicted, na.rm = T))
  return(pred.ranges)
}

auc_func<- function(test.data, predicted) {
  library(ROCR)
  temp<- dplyr::select(test.data, one_of(c("SED.TYPE", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS")))
  test.data<- data.frame(na.omit(temp))
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
    return(performance(dat, measure = "auc")@y.values[[1]])
  }
}

rmse_func<- function(test.data, predicted, response) {
  if(response == "Presence"){
    test.data<- dplyr::select(test.data, one_of(c("SED.TYPE", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    if(all(test.data$PRESENCE.BIOMASS == 0)){
      return(NA)
    } else {
      dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
      return(performance(dat, measure = "rmse")@y.values[[1]])
    }
  }
  
  if(response == "Biomass"){
    test.data<- dplyr::select(test.data, one_of(c("SED.TYPE", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "BIOMASS")
    if(all(test.data$BIOMASS == 0)){
      return(NA)
    } else {
      return(nrmse(sim = as.numeric(predicted), obs = test.data$BIOMASS))
    }
  }
}

calib_stat_func<- function(test.data, predicted){
  test.data<- dplyr::select(test.data, one_of(c("SED.TYPE", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
  test.data<- data.frame(na.omit(test.data))
  col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    calib.stat<- round(calibration(x = test.data[,col.ind], p = predicted)@statistic, 3)
    return(calib.stat)
  }
}

calib_plot_func<- function(test.data, predicted){
  test.data<- dplyr::select(test.data, one_of(c("SED.TYPE", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
  test.data<- data.frame(na.omit(test.data))
  col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    calib.all<- calibration(x = test.data[,col.ind], p = predicted)
    calib.df<- calib.all@calibration
    calib.plot<- ggplot() +
      geom_point(data = calib.df, aes(x = pedicted_pobability_bin, y = observed_poportion)) +
      xlim(c(0, 1)) + 
      ylim(c(0, 1)) +
      ylab("Proportion of Observed Occurrences") +
      xlab("Predicted Probability of Occurrence") +
      geom_abline(intercept = 0, linetype = "dashed")
    return(calib.plot)
  }
}

dat.full<- dat.full %>%
  left_join(., dat.test, by = c("COMNAME", "SEASON"))

mod.results<- data.frame(dat.full[,c(1:2)])
mod.results$DevExp.P<- rep(NA, nrow(mod.results))
mod.results$DevExp.B<- rep(NA, nrow(mod.results))
mod.results$Temp.DevExp.P<- rep(NA, nrow(mod.results))
mod.results$Temp.DevExp.B<- rep(NA, nrow(mod.results))
mod.results$AUC.SDM<- rep(NA, nrow(mod.results))
mod.results$RMSE.SDM.P<- rep(NA, nrow(mod.results))
mod.results$RMSE.SDM.B<- rep(NA, nrow(mod.results))
mod.results$Calib.SDM<- rep(NA, nrow(mod.results))
mod.results$Pred.Min.SDM<- rep(NA, nrow(mod.results))
mod.results$Pred.Max.SDM<- rep(NA, nrow(mod.results))
mod.results$Pred.BaseRate.SDM<- rep(NA, nrow(mod.results))
mod.results$RMSE.NEVA.B<- rep(NA, nrow(mod.results))

for(i in 1:nrow(dat.train)){
  
  # Get the data and fit the model
  dat.use<- dat.full[i,]
  season.use<- as.character(dat.use$SEASON[[1]])
  
  if(is.null(dat.use$Dir.Data[[1]])){
    print("No Vuln Assessment")
    next
  }
  
  dat.use<- dat.use %>%
    mutate(., "mod.fitted.p" = pmap(list(df = TRAIN.DATA, response = list("Presence")), possibly(gam_fit_full_func, NA)),
           "mod.fitted.b" = pmap(list(df = TRAIN.DATA, response = list("Biomass")), possibly(gam_fit_full_func, NA)))
  gam.mod0.p<- dat.use$mod.fitted.p[[1]]
  
  # No SDM data..
  if(is.infinite(summary(gam.mod0.p)$dev.expl) | summary(gam.mod0.p)$dev.expl > 0.97 | is.na(dat.use$mod.fitted.b)){
    print("Bad SDM fit")
    next()
  }
  
  # Else keep going
  saveRDS(dat.use$mod.fitted.p[[1]], file = paste(out.dir.stem, "gamfitpres", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
  saveRDS(dat.use$mod.fitted.b[[1]], file = paste(out.dir.stem, "gamfitbio", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
  
  # Store deviance explained
  mod.results$DevExp.P[i]<- round(summary(dat.use$mod.fitted.p[[1]])$dev.expl, 3)
  mod.results$DevExp.B[i]<- round(summary(dat.use$mod.fitted.b[[1]])$dev.expl, 3)
  
  # Compare with reduced model...
  gam.p.red<- gam(PRESENCE.BIOMASS ~ SED.TYPE + s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)
  gam.b.red<- gam(BIOMASS.MOD ~ SED.TYPE + s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = gaussian, select = TRUE)
  dev.p.temp<- (summary(dat.use$mod.fitted.p[[1]])$dev.expl - summary(gam.p.red)$dev.expl)/summary(dat.use$mod.fitted.p[[1]])$dev.expl
  dev.b.temp<- (summary(dat.use$mod.fitted.b[[1]])$dev.expl - summary(gam.b.red)$dev.expl)/summary(dat.use$mod.fitted.b[[1]])$dev.expl
  mod.results$Temp.DevExp.P[i]<- round(dev.p.temp, 3)
  mod.results$Temp.DevExp.B[i]<- round(dev.b.temp, 3)
  
  # Predictions and residuals
  dat.use<- dat.use %>%
    mutate(., "Predicted.SDM.P" = pmap(list(mod.fitted.p = mod.fitted.p, mod.fitted.b = mod.fitted.b, percentile = "SEASONALMU.OISST.Scale", test.data = TEST.DATA, response = "Presence"), possibly(predict_func, NA)),
           "Predicted.SDM.B" = pmap(list(mod.fitted = mod.fitted.p, mod.fitted.b = mod.fitted.b, percentile = "SEASONALMU.OISST.Scale", test.data = TEST.DATA, response = "Biomass"), possibly(predict_func, NA)),
           "Residuals.SDM.P" = pmap(list(test.data = TEST.DATA, response = list("Presence"), predicted = Predicted.SDM.P, type = "response"), possibly(resids_map_func, NA)),
           "Residuals.SDM.B" = pmap(list(test.data = TEST.DATA, response = list("Biomass"), predicted = Predicted.SDM.B, type = "response"), possibly(resids_map_func, NA)),
           "Pred.Ranges.SDM" = purrr::map(Predicted.SDM.P, possibly(pred_ranges_func, NA)),
           "AUC.SDM" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P), possibly(auc_func, NA)),
           "RMSE.SDM.P" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P, response = "Presence"), possibly(rmse_func, NA)),
           "RMSE.SDM.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.B, response = "Biomass"), possibly(rmse_func, NA)),
           "Calib.Stat.SDM" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P), possibly(calib_stat_func, NA)),
           "Calib.Plot.SDM" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P), possibly(calib_plot_func, NA)))
  
  mod.results$AUC.SDM[i]<- round(as.numeric(dat.use$AUC.SDM[[1]]), 3)
  mod.results$RMSE.SDM.P[i]<- round(as.numeric(dat.use$RMSE.SDM.P[[1]]), 3)
  mod.results$RMSE.SDM.B[i]<- round(as.numeric(dat.use$RMSE.SDM.B[[1]]), 3)
  mod.results$Calib.SDM[i]<- round(as.numeric(dat.use$Calib.Stat.SDM[[1]]), 3)
  mod.results$Pred.Min.SDM[i]<- round(as.numeric(dat.use$Pred.Ranges.SDM[[1]]$Min.Pred), 3)
  mod.results$Pred.Max.SDM[i]<- round(as.numeric(dat.use$Pred.Ranges.SDM[[1]]$Max.Pred), 3)
  mod.results$Pred.BaseRate.SDM[i]<- round(as.numeric(dat.use$Pred.Ranges.SDM[[1]]$Mean.Pred), 3)
  
  # A few plots --
  # Spatial residual maps
  resids.plot<- dat.use$Residuals.SDM.P[[1]]
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_PresResidsSDMMap.jpg", sep = ""), resids.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  resids.plot<- dat.use$Residuals.SDM.B[[1]]
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_BioResidsSDMMap.jpg", sep = ""), resids.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Validation plots
  calib.plot<- dat.use$Calib.Plot.SDM[[1]]
  if(any(is.na(calib.plot))){
    next
  }
  
  calib.plot<- calib.plot +
    annotate("text", x = 0.85, y = 0.15, label = paste("Pred ranges = c(", round(mod.results$Pred.Min.SDM[i], 2), ", ", round(mod.results$Pred.Max.SDM[i], 2), ")", sep = "")) +
    annotate("text", x = 0.85, y = 0.1, label = paste("Base rate = ", round(mod.results$Pred.BaseRate.SDM[i], 2))) +
    annotate("text", x = 0.85, y = 0.05, label = paste("AUC = ", round(mod.results$AUC.SDM[i], 2))) +
    annotate("text", x = 0.85, y = 0.0, label = paste("Calib = ", round(mod.results$Calib.SDM[i], 2)))
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_CalibSDMPlot.jpg", sep = ""), calib.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Draws of candidate models from the posterior for biomass component
  coef.mod0.b<- coef(dat.use$mod.fitted.b[[1]])
  vc.mod0.b<- vcov(dat.use$mod.fitted.b[[1]])
  cand.params.all.b<- rmvnorm(mod.sims, mean = coef.mod0.b, sigma = vc.mod0.b) 
  
  # Qualitative Data
  dir.dat<- dat.use$Dir.Data[[1]]
  sens.dat<- dat.use$Sens.Data[[1]]
  exp.dat<- dat.use$Exp.Data[[1]]
  
  # Adjustments
  sens.wt<- 5/5
  exp.wt<- 5/4
  dir.wt<- 5/3
  
  # Fix any of these?
  if(!is.null(fix.params)){
    if(fix.params == "All"){
      # Modify to match all column names of cand.params except intercept
      fix.params.use<- names(coef.mod0)[-which(names(coef.mod0) %in% "(Intercept)")]
      
      # Cand.params.all index...
      cand.params.all.ind<- which(colnames(cand.params.all) %in% fix.params.use, arr.ind = T)
      
      # Coef index
      coef.mod0.fix<- coef.mod0[names(coef.mod0) %in% fix.params.use]
      coef.mod0.fix.mat<- matrix(coef.mod0.fix, nrow = 1, ncol = length(coef.mod0.fix), byrow = T)
      
      # Move over values
      cand.params.all[,cand.params.all.ind]<- coef.mod0.fix.mat[rep(1:nrow(coef.mod0.fix.mat), times = nrow(cand.params.all)),]
      
      # Randomly generate intercept values...
      cand.params.all[,1]<- runif(mod.sims, -3, 3)
    } else {
      # Modify to match column names of cand.params
      fix.params.use<- paste(paste("s(", fix.params, ")", sep = ""), ".", seq(from = 1, to = 9, by = 1), sep = "")
      
      # Cand.params.all index...
      cand.params.all.ind<- which(colnames(cand.params.all) %in% fix.params.use, arr.ind = T)
      
      # Coef index
      coef.mod0.fix<- coef.mod0[names(coef.mod0) %in% fix.params.use]
      coef.mod0.fix.mat<- matrix(coef.mod0.fix, nrow = 1, ncol = length(coef.mod0.fix), byrow = T)
      
      # Move over values
      cand.params.all[,cand.params.all.ind]<- coef.mod0.fix.mat[rep(1:nrow(coef.mod0.fix.mat), times = nrow(cand.params.all)),]
    }
  }
  
  out.likes<- array(dim = c(mod.sims, 3))
  out.cands.b<- array(dim = c(mod.sims, ncol(cand.params.all.b)))
  
  for(k in 1:nrow(cand.params.all.b)){
    cand.params.b<- cand.params.all.b[k,]
    out.cands.b[k,]<- cand.params.b
    
    if(scenarios.type == "Dir") {
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0), nevaS = NULL, nevaE = NULL, dir.brks = dir.brks, fix.params = fix.params)
    } 
    
    if(scenarios.type == "Vuln"){
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = NULL, nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0), nevaE = round(c(exp.dat$Low, exp.dat$Moderate, exp.dat$High, exp.dat$Very.High)*exp.wt, 0), dir.brks = dir.brks, fix.params = fix.params)
    }
    
    if(scenarios.type == "Both"){
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0), nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0), nevaE = round(c(exp.dat$Low, exp.dat$Moderate, exp.dat$High, exp.dat$Very.High)*exp.wt, 0), dir.brks = dir.brks, fix.params = fix.params)
    }
    
    if(scenarios.type == "Dir.Sens"){
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0), nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0), nevaE = NULL, dir.brks = dir.brks, fix.params = fix.params)
    }
  }
  
  # Bring together results in a list
  out.likes[is.infinite(out.likes)]<- NA
  mcmc.results<- list(out.likes, out.cands.b)
  
  # Processing and results plots
  likes.df<- data.frame(mcmc.results[[1]])
  likes.nsamps<- nrow(likes.df)
  
  posts.df.b<- data.frame(mcmc.results[[2]])
  posts.nsamps<- nrow(posts.df.b)
  
  # Likes dataframe and plots of mixing
  names(likes.df)<- c("Likelihood (NEVA|SDM)", "Prior P(SDM.B)", "Posterior")
  likes.df$Iteration<- rep(seq(from = 1, to = likes.nsamps, by = 1))
  likes.df<- likes.df %>%
    gather(., Sample, Value, -Iteration)
  likes.df$Sample<- factor(likes.df$Sample, levels = c("Likelihood (NEVA|SDM)", "Prior P(SDM.B)", "Posterior"))
  
  out.plot<- ggplot(likes.df) +
    geom_line(aes(x = Iteration, y = Value), alpha = 0.25) +
    #scale_color_manual(values = c('#4daf4a', '#377eb8', '#e41a1c')) + 
    ylab("Log Likelihood") + 
    facet_wrap(~Sample, scales = "free") + 
    theme_bw()
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_LikePriorPost.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Predictions from best model 
  like.posts.b<- likes.df[likes.df$Sample == "Prior P(SDM.B)",]
  posts.samps.b<- posts.df.b
  
  # Get maximimum value and extract row from posts.samps
  mod.ind.b<- like.posts.b$Iteration[which.max(like.posts.b$Value)]
  best.fit.b<- posts.samps.b[mod.ind.b,]
  best.fit.mat.b<- matrix(as.numeric(best.fit.b), nrow = 1, ncol = length(best.fit.b), byrow = T, dimnames = list(NULL, names(coef(dat.use$mod.fitted.b))))
  
  # Make predictions with these values
  dat.use<- dat.use %>%
    mutate(., "Predicted.NEVA.B" = pmap(list(mod.fitted.p = mod.fitted.p, mod.fitted.b = mod.fitted.b, test.data = TEST.DATA, neva.best.fit.b = list(best.fit.mat.b)), possibly(neva_predict_bio_func, NA)),
           "Residuals.NEVA.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.NEVA.B, response = "Biomass", type = "response"), possibly(resids_map_func, NA)),
           "RMSE.NEVA.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.NEVA.B, response = "Biomass"), possibly(rmse_func, NA)))
  
  # Store key results
  mod.results$RMSE.NEVA.B[i]<- round(as.numeric(dat.use$RMSE.NEVA.B[[1]]), 3)
  
  # A few plots --
  # Spatial residual maps
  resids.plot<- dat.use$Residuals.NEVA.B[[1]]
  ggplot2::ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_BioResidsNEVAMap.jpg", sep = ""), resids.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Prediction Maps
  source("~/GitHub/COCA/Code/plot_func_06132018.R")
  maps<- plot_func(response = "Presence", mod.fitted.p = dat.use$mod.fitted.p[[1]], mod.fitted.b = dat.use$mod.fitted.b[[1]], season = season.use, base.preds = base.preds, fut.preds = fut.preds, like.posts = like.posts.b, posts.samps = posts.df.b)
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_PresenceMaps.jpg", sep = ""), maps, width = 10, height = 11, dpi = 125, units = "in")
  rm(maps)
  maps<- plot_func(response = "Biomass", mod.fitted.p = dat.use$mod.fitted.p[[1]], mod.fitted.b = dat.use$mod.fitted.b[[1]], season = season.use, base.preds = base.preds, fut.preds = fut.preds, like.posts = like.posts.b, posts.samps = posts.df.b)
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_BiomassMaps.jpg", sep = ""), maps, width = 10, height = 11, dpi = 125, units = "in")
  rm(maps)
  
  # All possible prediction curves, the original smooth, and then the one we have selected...
  # Some prep
  ilink<- family(dat.use$mod.fitted.b[[1]])$linkinv
  pred.dat.use<- pred.dat[pred.dat$SEASON == season.use, ]
  rescaled.dat.use<- rescaled.dat[rescaled.dat$Season == season.use, ]
  
  # Predictor matrix
  pred.mat<- predict(dat.use$mod.fitted.b[[1]], newdata = pred.dat.use, type = "lpmatrix")
  
  # Original predictions
  pred0<- exp(ilink(pred.mat %*% coef(dat.use$mod.fitted.b[[1]])))
  
  # Best prediction
  like.posts = likes.df[likes.df$Sample == "Posterior",]
  mod.ind<- like.posts$Iteration[which.max(like.posts$Value)]
  pred.best.vec<- as.numeric(posts.df.b[mod.ind,])
  names(pred.best.vec)<- names(coef(dat.use$mod.fitted.b[[1]]))
  pred.best<- exp(ilink(pred.mat %*% pred.best.vec))
  
  # All predictions
  pred.all<- apply(cand.params.all.b, 1, function(x) exp(ilink(pred.mat %*% x)))
  
  # Plotting
  # SST
  want<- 1:500
  sst.all<- data.frame(pred.all[want,])
  colnames(sst.all)<- paste("Mod.", seq(from = 1, to = ncol(sst.all)), sep = "")
  sst.all<- sst.all %>%
    gather(., "Model", "Pred")
  sst.all$Value<- rescaled.dat.use$SST
  sst.all$Parameter<- rep("SST", nrow(sst.all))
  
  sst.base<- ggplot(sst.all, aes(x = Value, y = Pred, group = Model)) + 
    geom_line(show.legend = F, color = "#bdbdbd") 
  
  sst.dat<- data.frame("Parameter" = rep("SST", 1000), "Value" = rep(rescaled.dat.use$SST, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  sst.out<- sst.base +
    geom_line(data = sst.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("SST") +
    theme_bw() 
  
  # Depth
  want<- 501:1000
  depth.all<- data.frame(pred.all[want,])
  colnames(depth.all)<- paste("Mod.", seq(from = 1, to = ncol(depth.all)), sep = "")
  depth.all<- depth.all %>%
    gather(., "Model", "Pred")
  depth.all$Value<- rescaled.dat.use$Depth
  depth.all$Parameter<- rep("Depth", nrow(depth.all))
  
  depth.base<- ggplot(depth.all, aes(x = Value, y = Pred, group = Model)) + 
    geom_line(show.legend = F, color = "#bdbdbd")
  
  depth.dat<- data.frame("Parameter" = rep("Depth", 1000), "Value" = rep(rescaled.dat.use$Depth, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  depth.out<- depth.base +
    geom_line(data = depth.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("Depth") +
    theme_bw() 
  
  out<- plot_grid(sst.out + theme(legend.position="none"), depth.out + theme(legend.position="none"), nrow = 1, ncol = 2, align = "hv", scale = 1)
  legend<- get_legend(sst.out)
  out<- plot_grid(out, legend, rel_widths = c(3, 0.5))
  ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_PredCurv.jpg", sep = ""), out, width = 11, height = 8, dpi = 125, units = "in")
  
  # Update
  print(paste(tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), " is done!", sep = ""))
  
  # Save MCMC results
  saveRDS(mcmc.results, file = paste(out.dir.stem, "mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
}

write.csv(mod.results, file = paste(out.dir.stem, "clam.mod.results.csv", sep = ""))

#####
## SDM predictions
#####
res.files<- list.files(out.dir.stem, paste("mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), sep = ""))

for(i in seq_along(res.files)){
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  spp.season.match<- paste(spp, season, sep = ".")
  
  temp.dev.exp<- mod.results$Temp.DevExp.B[paste(mod.results$COMNAME, mod.results$SEASON, sep = ".") == spp.season.match]
  
  # Model fit -- presence and biomass SDM
  mod.fitted.p<- readRDS(paste(out.dir.stem, gsub("mcmc_", "gamfitpres", res.files[i]), sep = ""))
  mod.fitted.b<- readRDS(paste(out.dir.stem, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(mod.fitted.b)$linkinv
  gam.coef<- names(coef(mod.fitted.p))
  
  # Model fit -- NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  mod.ind.b<- likes.temp$Iteration[which.max(likes.temp$Posterior)]
  
  posts.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[2]])
  best.fit<- posts.temp[mod.ind.b,]
  best.fit.mat<- matrix(as.numeric(best.fit), nrow = 1, ncol = length(best.fit), byrow = T, dimnames = list(NULL, gam.coef))
  
  # Make predictions
  # SDM
  sdm.base<- base.preds$Data[[match(season, base.preds$SEASON)]]
  sdm.base<- sdm.base %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
  sdm.base.p<- predict.gam(mod.fitted.p, newdata = sdm.base, type = "response")
  sdm.base.b<- as.numeric(sdm.base.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = sdm.base, type = "response")))
  
  # Mean climate model SST
  newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  newdat.mu<- newdat.mu %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
  sdm.fut.mu.p<- predict.gam(mod.fitted.p, newdata = newdat.mu, type = "response")
  sdm.fut.mu.b<- as.numeric(sdm.fut.mu.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.mu, type = "response")))
  
  # pct05 climate model
  newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.05<- newdat.05 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
  names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
  sdm.fut.pct05.p<- predict.gam(mod.fitted.p, newdata = newdat.05, type = "response")
  sdm.fut.pct05.b<- as.numeric(sdm.fut.pct05.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.05, type = "response")))
  
  # pct95 climate model
  newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.95<- newdat.95 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
  names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
  sdm.fut.pct95.p<- predict.gam(mod.fitted.p, newdata = newdat.95, type = "response")
  sdm.fut.pct95.b<- as.numeric(sdm.fut.pct95.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.95, type = "response")))
  
  # Presence only dataframes
  sdm.map.base.p<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base.p)
  sdm.map.fut.mu.p<- data.frame("x" = newdat.mu$x, "y" = newdat.mu$y, "pred" = sdm.fut.mu.p)
  sdm.map.fut.pct05.p<- data.frame("x" = newdat.05$x, "y" = newdat.05$y, "pred.05" = sdm.fut.pct05.p)
  sdm.map.fut.pct95.p<- data.frame("x" = newdat.95$x, "y" = newdat.95$y, "pred.95" = sdm.fut.pct95.p)
  sdm.map.base.b<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base.b)
  sdm.map.fut.mu.b<- data.frame("x" = newdat.mu$x, "y" = newdat.mu$y, "pred" = sdm.fut.mu.b)
  sdm.map.fut.pct05.b<- data.frame("x" = newdat.05$x, "y" = newdat.05$y, "pred.05" = sdm.fut.pct05.b)
  sdm.map.fut.pct95.b<- data.frame("x" = newdat.95$x, "y" = newdat.95$y, "pred.95" = sdm.fut.pct95.b)
  
  # Differences
  sdm.diff.p<- data.frame("x" = sdm.map.fut.mu.p$x, "y" = sdm.map.fut.mu.p$y, "pred" = sdm.map.fut.mu.p$pred - sdm.map.base.p$pred)
  sdm.lwr.diff.p<- data.frame("x" = sdm.map.fut.pct05.p$x, "y" = sdm.map.fut.pct05.p$y, "pred" = sdm.map.fut.pct05.p$pred - sdm.map.base.p$pred)
  sdm.upr.diff.p<- data.frame("x" = sdm.map.fut.pct95.p$x, "y" = sdm.map.fut.pct95.p$y, "pred" = sdm.map.fut.pct95.p$pred - sdm.map.base.p$pred)
  names(sdm.map.base.p)[3]<- "Baseline.sdm.p"
  names(sdm.map.fut.mu.p)[3]<- "Future_mean.sdm.p"
  names(sdm.map.fut.pct05.p)[3]<- "Future_cold.sdm.p"
  names(sdm.map.fut.pct95.p)[3]<- "Future_warm.sdm.p"
  names(sdm.diff.p)[3]<- "Future_mean_diff.sdm.p"
  names(sdm.lwr.diff.p)[3]<- "Future_cold_diff.sdm.p"
  names(sdm.upr.diff.p)[3]<- "Future_warm_diff.sdm.p"
  
  sdm.diff.b<- data.frame("x" = sdm.map.fut.mu.b$x, "y" = sdm.map.fut.mu.b$y, "pred" = sdm.map.fut.mu.b$pred - sdm.map.base.b$pred)
  sdm.lwr.diff.b<- data.frame("x" = sdm.map.fut.pct05.b$x, "y" = sdm.map.fut.pct05.b$y, "pred" = sdm.map.fut.pct05.b$pred - sdm.map.base.b$pred)
  sdm.upr.diff.b<- data.frame("x" = sdm.map.fut.pct95.b$x, "y" = sdm.map.fut.pct95.b$y, "pred" = sdm.map.fut.pct95.b$pred - sdm.map.base.b$pred)
  names(sdm.map.base.b)[3]<- "Baseline.sdm.b"
  names(sdm.map.fut.mu.b)[3]<- "Future_mean.sdm.b"
  names(sdm.map.fut.pct05.b)[3]<- "Future_cold.sdm.b"
  names(sdm.map.fut.pct95.b)[3]<- "Future_warm.sdm.b"
  names(sdm.diff.b)[3]<- "Future_mean_diff.sdm.b"
  names(sdm.lwr.diff.b)[3]<- "Future_cold_diff.sdm.b"
  names(sdm.upr.diff.b)[3]<- "Future_warm_diff.sdm.b"
  
  # Getting combined biomass predictions
  # Make predictions with these values
  lpmat.base<- predict.gam(mod.fitted.b, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "lpmatrix")
  combo.map.base<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base.p * exp(ilink(as.numeric(lpmat.base %*% t(best.fit.mat)))))
  
  # Mean climate model SST
  newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  newdat.mu<- newdat.mu %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
  lpmat.fut.mu<- predict.gam(mod.fitted.b, newdata = newdat.mu, type = "lpmatrix")
  combo.map.fut.mu<- data.frame("x" = newdat.mu$x, "y" = newdat.mu$y, "pred" = sdm.fut.mu.p * exp(ilink(as.numeric(lpmat.fut.mu %*% t(best.fit.mat)))))
  
  # pct05 climate model
  newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.05<- newdat.05 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
  names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
  lpmat.fut.pct05<- predict.gam(mod.fitted.b, newdata = newdat.05, type = "lpmatrix")
  combo.map.fut.pct05<- data.frame("x" = newdat.05$x, "y" = newdat.05$y, "pred.05" = sdm.fut.pct05.p * exp(ilink(as.numeric(lpmat.fut.pct05 %*% t(best.fit.mat)))))
  
  # pct95 climate model
  newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.95<- newdat.95 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
  names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
  lpmat.fut.pct95<- predict.gam(mod.fitted.b, newdata = newdat.95, type = "lpmatrix")
  combo.map.fut.pct95<- data.frame("x" = newdat.95$x, "y" = newdat.95$y, "pred.95" = sdm.fut.pct95.p * exp(ilink(as.numeric(lpmat.fut.pct95 %*% t(best.fit.mat)))))
  
  # Differences
  combo.diff<- data.frame("x" = combo.map.fut.mu$x, "y" = combo.map.fut.mu$y, "pred" = combo.map.fut.mu$pred - combo.map.base$pred)
  combo.lwr.diff<- data.frame("x" = combo.map.fut.pct05$x, "y" = combo.map.fut.pct05$y, "pred" = combo.map.fut.pct05$pred - combo.map.base$pred)
  combo.upr.diff<- data.frame("x" = combo.map.fut.pct95$x, "y" = combo.map.fut.pct95$y, "pred" = combo.map.fut.pct95$pred - combo.map.base$pred)
  
  # Calculate fish availability for the different datasets: combo.map.base, combo.map.fut.mu, combo.map.fut.pct05, combo.map.fut.pct95, combo.diff, combo.lwr.diff, combo.upr.diff
  names(combo.map.base)[3]<- "Baseline.combo.b"
  names(combo.map.fut.mu)[3]<- "Future_mean.combo.b"
  names(combo.map.fut.pct05)[3]<- "Future_cold.combo.b"
  names(combo.map.fut.pct95)[3]<- "Future_warm.combo.b"
  names(combo.diff)[3]<- "Future_mean_diff.combo.b"
  names(combo.lwr.diff)[3]<- "Future_cold_diff.combo.b"
  names(combo.upr.diff)[3]<- "Future_warm_diff.combo.b"
  
  projections.dat<- sdm.map.base.p %>%
    left_join(., sdm.map.fut.mu.p, by = c("x", "y")) %>%
    left_join(., sdm.map.fut.pct05.p, by = c("x", "y")) %>%
    left_join(., sdm.map.fut.pct95.p, by = c("x", "y")) %>%
    left_join(., sdm.diff.p, by = c("x", "y")) %>%
    left_join(., sdm.lwr.diff.p, by = c("x", "y")) %>%
    left_join(., sdm.upr.diff.p, by = c("x", "y")) %>%
    left_join(., sdm.map.base.b, by = c("x", "y")) %>%
    left_join(., sdm.map.fut.mu.b, by = c("x", "y")) %>%
    left_join(., sdm.map.fut.pct05.b, by = c("x", "y")) %>%
    left_join(., sdm.map.fut.pct95.b, by = c("x", "y")) %>%
    left_join(., sdm.diff.b, by = c("x", "y")) %>%
    left_join(., sdm.lwr.diff.b, by = c("x", "y")) %>%
    left_join(., sdm.upr.diff.b, by = c("x", "y")) %>%
    left_join(., combo.map.base, by = c("x", "y")) %>%
    left_join(., combo.map.fut.mu, by = c("x", "y")) %>%
    left_join(., combo.map.fut.pct05, by = c("x", "y")) %>%
    left_join(., combo.map.fut.pct95, by = c("x", "y")) %>%
    left_join(., combo.diff, by = c("x", "y")) %>%
    left_join(., combo.lwr.diff, by = c("x", "y")) %>%
    left_join(., combo.upr.diff, by = c("x", "y")) 
  projections.dat$COMNAME<- spp
  projections.dat$SEASON<- season
  projections.dat<- projections.dat %>%
    gather(., "Proj.Class", "Projection", -x, -y, -COMNAME, -SEASON) %>%
    group_by(., COMNAME, SEASON, Proj.Class) %>%
    nest(.key = "Projections")
  
  if(i == 1){
    result<- projections.dat
    print(paste(spp, season, "is done!", sep = " "))
  } else {
    result<- bind_rows(result, projections.dat)
    print(paste(spp, season, "is done!", sep = " "))
  }
}

saveRDS(result, file = paste(out.dir.stem, "atlantic surfclam", "SDMPredictions.rds", sep = ""))

######
## SDM vs. Expert expectations: Sensitivity = magnitude of change, Directional effect = net change
######
jon.df<- read_csv("~/GitHub/COCA/Data/Jon_QualitativeResults.csv")
result<- read_rds(paste(out.dir.stem, "atlantic surfclam", "SDMPredictions.rds", sep = ""))
result<- result %>%
  left_join(., jon.df, by = "COMNAME")

# Interested in SDM presence, SDM biomass
comps<- c("Future_mean_diff.sdm.p", "Future_mean_diff.sdm.b")

means_func<- function(df){
  mean(df$Projection, na.rm = T)
}

abs_func<- function(means){
  abs(means)
}

result.comps<- result %>%
  dplyr::filter(., Proj.Class %in% comps)

for(i in seq_along(comps)){
  plot.df<- result.comps %>%
    dplyr::filter(., Proj.Class == comps[i]) %>%
    mutate(., "Mean.Diff" = as.numeric(map(Projections, means_func)))
  
  plot.df$DIRECTIONAL.EFFECT<- factor(plot.df$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"))
  plot.df$DIRECTIONAL.EFFECT.CERTAINTY<- factor(plot.df$DIRECTIONAL.EFFECT.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))
  plot.df$SENSITIVITY<- factor(plot.df$SENSITIVITY, levels = c("Low", "Moderate", "High", "Very high"))
  plot.df$SENSITIVITY.CERTAINTY<- factor(plot.df$SENSITIVITY.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))
  
  # Overall certainty....
  # Join with overall uncertainty...
  certainty.levels<- data.frame("SENSITIVITY.CERTAINTY" = c("Low", "Moderate", "High", "Very high"), "DIRECTIONAL.EFFECT.CERTAINTY" = c("Low", "Moderate", "High", "Very high"), "Score" = c(1, 2, 3, 4))
  
  # Jon's data
  plot.df<- plot.df %>%
    mutate(., "Certainty" = paste(SENSITIVITY.CERTAINTY, "_", DIRECTIONAL.EFFECT.CERTAINTY, sep = ""))
  
  plot.df<- plot.df %>%
    left_join(., dplyr::select(certainty.levels, c(SENSITIVITY.CERTAINTY, Score)), by = "SENSITIVITY.CERTAINTY")
  names(plot.df)[13]<- "SENSITIVITY.Score"
  
  plot.df<- plot.df %>%
    left_join(., dplyr::select(certainty.levels, c(DIRECTIONAL.EFFECT.CERTAINTY, Score)), by = "DIRECTIONAL.EFFECT.CERTAINTY")
  names(plot.df)[14]<- "DirectionalEffect.Score"
  
  plot.df$Certainty.Score<- plot.df$SENSITIVITY.Score + plot.df$DirectionalEffect.Score
  plot.df$Certainty.Score<- factor(plot.df$Certainty.Score, levels = c(2, 3, 4, 5, 6, 7, 8), labels = c("Low", "Low_Moderate", "Moderate", "Moderate_High", "High", "High_Very high", "Very high"))
  
  plot.df$DIRECTIONAL.EFFECT<- factor(plot.df$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"))
  plot.df$DIRECTIONAL.EFFECT.CERTAINTY<- factor(plot.df$DIRECTIONAL.EFFECT.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))
  plot.df$SENSITIVITY<- factor(plot.df$SENSITIVITY, levels = c("Low", "Moderate", "High", "Very high"))
  plot.df$SENSITIVITY.CERTAINTY<- factor(plot.df$SENSITIVITY.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))
  
  plot.df<- plot.df %>%
    mutate(., "Future_mean_diff_abs" = as.numeric(map(Mean.Diff, abs_func)))
  
  plot.df.avg<- plot.df %>%
    group_by(., SEASON, SENSITIVITY.CERTAINTY, SENSITIVITY) %>%
    summarise_at(., "Future_mean_diff_abs", c("mean", "sd"), na.rm = T) %>%
    drop_na(SENSITIVITY)
  
  plot.out.sens<- ggplot(data = plot.df.avg, aes(x = SENSITIVITY, y = mean)) +
    geom_point() + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    ylab("Average absolute change") +
    xlab("Climate sensitivity") +
    facet_wrap(~SEASON+SENSITIVITY.CERTAINTY, nrow = 2, scales = "free_y") +
    geom_smooth(data = plot.df.avg, aes(x = as.numeric(SENSITIVITY), y = mean, color = SENSITIVITY.CERTAINTY), method = "lm", se = FALSE) +
    scale_color_manual(name = "Sensitivity certainty", values = c('#b2e2e2','#66c2a4','#2ca25f','#006d2c')) +
    theme_bw() 
  plot.title<- gsub("Future_mean_diff.", "", comps[i])
  ggsave(paste(out.dir.stem, "Sensitivity_vs_Absolute", plot.title, "Change.jpg", sep = ""), width = 11, height = 8, units = "in")
  
  # Now with directional effect
  plot.df.avg<- plot.df %>%
    group_by(., SEASON, DIRECTIONAL.EFFECT.CERTAINTY, DIRECTIONAL.EFFECT) %>%
    summarise_at(., "Mean.Diff", c("mean", "sd"), na.rm = T) %>%
    drop_na(DIRECTIONAL.EFFECT)
  
  plot.out.dir<- ggplot(data = plot.df.avg, aes(x = DIRECTIONAL.EFFECT, y = mean)) +
    geom_point() + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    ylab("Average change") +
    xlab("Directional effect") +
    facet_wrap(~SEASON+DIRECTIONAL.EFFECT.CERTAINTY, nrow = 2, scales = "free_y") +
    geom_smooth(data = plot.df.avg, aes(x = as.numeric(DIRECTIONAL.EFFECT), y = mean, color = DIRECTIONAL.EFFECT.CERTAINTY), method = "lm", se = FALSE) +
    scale_color_manual(name = "Directional effect certainty", values = c('#b2e2e2','#66c2a4','#2ca25f','#006d2c')) +
    theme_bw() 
  ggsave(paste(out.dir.stem, "DirEffect_vs_", plot.title, "Change.jpg", sep = ""), width = 11, height = 8, units = "in")
  
}


######
## SDM+NEVA vs. SDM plotting of prediction intervals?
######
rmse_func<- function(predicted, test.data) {
  if(all(test.data$BIOMASS == 0)){
    return(NA)
  } else {
    return(nrmse(sim = as.numeric(as.data.frame(predicted)[,1]), obs = test.data$BIOMASS))
  }
}

# Going to need access to the likelihood, prior and posterior information. This should be all in the mcmc results file for each species...
res.files<- list.files(out.dir.stem, "mcmc_")

# Results
neva.vs.sdm<- data.frame("COMNAME" = rep(NA, length(res.files)), "SEASON" = rep(NA, length(res.files)), "NEVA.RMSE.mu" = rep(NA, length(res.files)), "NEVA.RMSE.sd" = rep(NA, length(res.files)), "SDM.RMSE.mu" = rep(NA, length(res.files)), "SDM.RMSE.sd" = rep(NA, length(res.files)))

# Model results -- should really only be thinking about this for "good" models? AUC greater than 0.7?
auc.res<- read_csv(paste(out.dir.stem, "mod.results.csv", sep = ""))

for(i in seq_along(res.files)){
  
  # Model fit -- presence and biomass
  gam.p.temp<- readRDS(paste(out.dir.stem, gsub("mcmc_", "gamfitpres", res.files[i]), sep = ""))
  gam.b.temp<- readRDS(paste(out.dir.stem, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(gam.b.temp)$linkinv
  gam.coef<- names(coef(gam.p.temp))
  
  # Ranks of candidate draws...
  # Get ranks for SDM only and for NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  likes.temp$SDM.Rank<- order(likes.temp$Prior)
  likes.temp$NEVA.Rank<- order(likes.temp$Posterior)
  
  # Need to make predictions from each of these iterations
  mods.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[2]])
  colnames(mods.temp)<- gam.coef
  
  # Predictions
  spp.match<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season.match<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  col.check<- paste(spp.match, season.match, sep = ".")
  
  # Keep going?
  auc.temp<- auc.res$AUC.SDM[paste(auc.res$COMNAME, auc.res$SEASON, sep = ".") == col.check] 
  print(auc.temp)
  
  if(auc.temp >= 0.7){
    test.data<- dat.full %>%
      filter(., COMNAME == spp.match & SEASON == season.match) %>%
      dplyr::select(TEST.DATA) %>%
      unnest() %>%
      data.frame()
    temp<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS")))
    test.data<- data.frame(na.omit(temp))
    
    # Presence component -- doesn't change
    pred.p<- predict.gam(gam.p.temp, newdata = test.data, type = "response")
    
    # Prediction storage
    preds.out<- data.frame("Pt" = seq(from = 1, nrow(test.data)), "Pred" = rep(NA, nrow(test.data)))
    
    for(j in 1:nrow(mods.temp)){
      fit.mat<- matrix(as.numeric(mods.temp[j,]), nrow = 1, ncol = length(mods.temp[j,]), byrow = T, dimnames = list(NULL, gam.coef))
      
      # Make predictions with these values
      lpmat.pred<- predict.gam(gam.b.temp, newdata = test.data, type = "lpmatrix")
      preds.out[,j+1]<- as.numeric(exp(ilink(as.numeric(lpmat.pred %*% t(fit.mat))))*pred.p)
      names(preds.out)[j+1]<- paste("Iteration.", j, sep = "")
    }
    
    # Get the SDM predictions...
    sdm.best<- likes.temp[likes.temp$SDM.Rank<= 900, ]
    sdm.best<- sdm.best[order(sdm.best$SDM.Rank), ]
    sdm.preds.best<- preds.out[,sdm.best$Iteration+1] %>%
      gather(., "Iteration", "Prediction") %>%
      separate(., Iteration, c("Remove", "Iteration")) %>%
      dplyr::select(., -Remove) %>%
      group_by(Iteration) %>%
      nest(., .key = "Pred")
    
    # Get the NEVA predictions
    neva.best<- likes.temp[likes.temp$NEVA.Rank<= 900, ]
    neva.best<- neva.best[order(neva.best$NEVA.Rank), ]
    neva.preds.best<- preds.out[,neva.best$Iteration+1] %>%
      gather(., "Iteration", "Prediction") %>%
      separate(., Iteration, c("Remove", "Iteration")) %>%
      dplyr::select(., -Remove) %>%
      group_by(Iteration) %>%
      nest(., .key = "Pred") 
    
    # Calculate the AUC and Calibration values
    test.data<- dat.full %>%
      filter(., COMNAME == spp.match & SEASON == season.match) %>%
      dplyr::select(TEST.DATA) %>%
      unnest() %>%
      data.frame()
    temp<- dplyr::select(test.data, one_of(c("BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale")))
    test.data<- data.frame(na.omit(temp))
    
    sdm.summary<- sdm.preds.best %>%
      mutate(., "RMSE" = map2(Pred, list(test.data), possibly(rmse_func, NA))) %>%
      dplyr::select(., Iteration, RMSE) %>%
      unnest() %>%
      summarize_at(., vars("RMSE"), funs(mean, sd), na.rm = T)
    
    neva.summary<- neva.preds.best %>%
      mutate(., "RMSE" = map2(Pred, list(test.data), possibly(rmse_func, NA))) %>%
      dplyr::select(., Iteration, RMSE) %>%
      unnest() %>%
      summarize_at(., vars("RMSE"), funs(mean, sd), na.rm = T)
    
    neva.vs.sdm$COMNAME[i]<- spp.match
    neva.vs.sdm$SEASON[i]<- season.match
    neva.vs.sdm$NEVA.RMSE.mu[i]<- neva.summary$mean
    neva.vs.sdm$NEVA.RMSE.sd[i]<- neva.summary$sd
    neva.vs.sdm$SDM.RMSE.mu[i]<- sdm.summary$mean
    neva.vs.sdm$SDM.RMSE.sd[i]<- sdm.summary$sd
  } else {
    print(paste(season, spp.match, "doesn't pass auc.check", sep = " "))
    neva.vs.sdm$COMNAME[i]<- spp.match
    neva.vs.sdm$SEASON[i]<- season.match
    neva.vs.sdm$NEVA.RMSE.mu[i]<- NA
    neva.vs.sdm$NEVA.RMSE.sd[i]<- NA
    neva.vs.sdm$SDM.RMSE.mu[i]<- NA
    neva.vs.sdm$SDM.RMSE.sd[i]<- NA
  }
}

# Did the AUC crop work?
auc.res$Col.Check<- paste(auc.res$COMNAME, auc.res$SEASON, sep = ".")
neva.vs.sdm$Col.Check<- paste(neva.vs.sdm$COMNAME, neva.vs.sdm$SEASON, sep = ".")
temp<- auc.res$AUC.SDM[auc.res$Col.Check %in% neva.vs.sdm$Col.Check]

rmse.cert.dat.mu<- neva.vs.sdm %>%
  dplyr::select(., COMNAME, SEASON, NEVA.RMSE.mu, SDM.RMSE.mu) %>%
  mutate(., "RMSE.mu.diff" = round(NEVA.RMSE.mu - SDM.RMSE.mu), 3) %>%
  dplyr::select(., COMNAME, SEASON, RMSE.mu.diff) 
rmse.cert.dat.sd<- neva.vs.sdm %>%
  dplyr::select(., COMNAME, SEASON, NEVA.RMSE.sd, SDM.RMSE.sd) %>%
  mutate(., "RMSE.sd.diff" = round(NEVA.RMSE.sd - SDM.RMSE.sd), 3) %>%
  dplyr::select(., COMNAME, SEASON, RMSE.sd.diff) 
rmse.cert.dat<- data.frame(rmse.cert.dat.mu, "RMSE.sd.diff" = rmse.cert.dat.sd$RMSE.sd.diff)

# Join with overall uncertainty...
certainty.levels<- data.frame("SENSITIVITY.CERTAINTY" = c("Low", "Moderate", "High", "Very high"), "DIRECTIONAL.EFFECT.CERTAINTY" = c("Low", "Moderate", "High", "Very high"), "Score" = c(1, 2, 3, 4))

# Jon's data
jon.df<- read_csv("~/GitHub/COCA/Data/Jon_QualitativeResults.csv")

rmse.cert.dat<- rmse.cert.dat %>%
  left_join(., jon.df, by = "COMNAME")

rmse.cert.dat<- rmse.cert.dat %>%
  mutate(., "Certainty" = paste(SENSITIVITY.CERTAINTY, "_", DIRECTIONAL.EFFECT.CERTAINTY, sep = ""))

rmse.cert.dat<- rmse.cert.dat %>%
  left_join(., dplyr::select(certainty.levels, c(SENSITIVITY.CERTAINTY, Score)), by = "SENSITIVITY.CERTAINTY")
names(rmse.cert.dat)[12]<- "SENSITIVITY.Score"

rmse.cert.dat<- rmse.cert.dat %>%
  left_join(., dplyr::select(certainty.levels, c(DIRECTIONAL.EFFECT.CERTAINTY, Score)), by = "DIRECTIONAL.EFFECT.CERTAINTY")
names(rmse.cert.dat)[13]<- "DirectionalEffect.Score"

rmse.cert.dat$Certainty.Score<- rmse.cert.dat$SENSITIVITY.Score + rmse.cert.dat$DirectionalEffect.Score
rmse.cert.dat$Certainty.Score<- factor(rmse.cert.dat$Certainty.Score, levels = c(2, 3, 4, 5, 6, 7, 8), labels = c("Low", "Low_Moderate", "Moderate", "Moderate_High", "High", "High_Very high", "Very high"))
rmse.cert.dat$RMSE.mu.diff.code<- ifelse(rmse.cert.dat$RMSE.mu.diff > 0, "Better", 
                                         ifelse(rmse.cert.dat$RMSE.mu.diff == 0, "No change", "Worse"))
rmse.cert.dat$RMSE.mu.diff.code<- factor(rmse.cert.dat$RMSE.mu.diff.code, levels = c("Better", "No change", "Worse"))
rmse.cert.dat$RMSE.sd.diff.code<- ifelse(rmse.cert.dat$RMSE.sd.diff > 0, "Better", 
                                         ifelse(rmse.cert.dat$RMSE.sd.diff == 0, "No change", "Worse"))
rmse.cert.dat$RMSE.sd.diff.code<- factor(rmse.cert.dat$RMSE.sd.diff.code, levels = c("Better", "No change", "Worse"))

rmse.cert.dat$SENSITIVITY.CERTAINTY<- factor(rmse.cert.dat$SENSITIVITY.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))
rmse.cert.dat$DIRECTIONAL.EFFECT.CERTAINTY<- factor(rmse.cert.dat$DIRECTIONAL.EFFECT.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))

rmse.cert.plot<- ggplot(na.omit(rmse.cert.dat), aes(x = SENSITIVITY.CERTAINTY, y = DIRECTIONAL.EFFECT.CERTAINTY, color = RMSE.mu.diff.code, label = str_to_title(COMNAME))) + 
  scale_color_manual(name = "Effect of including\nNEVA on RMSE", values = c("#4daf4a", "gray", "#e41a1c")) + 
  geom_text_repel(segment.color = NA, show.legend = F) +
  xlab("Sensitivity certainty") +
  ylab("Directional effect certainty") +
  facet_wrap(~SEASON, scales = "free", nrow = 2)
ggsave(paste(out.dir.stem, "rmse.plot.jpg", sep = ""), height = 8, width = 11, units = "in")


#####
## SDM+NEVA vs. SDM plotting of parameter intervals?
#####
# Going to need access to the likelihood, prior and posterior information. This should be all in the mcmc results file for each species...
res.files<- list.files(out.dir.stem, "mcmc_")

# Results
for(i in seq_along(res.files)){
  
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  
  # Model fit -- presence and biomass
  gam.b.temp<- readRDS(paste(out.dir.stem, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(gam.b.temp)$linkinv
  gam.coef<- names(coef(gam.b.temp))
  
  # Ranks of candidate draws...
  # Get ranks for SDM only and for NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  likes.temp$SDM.Rank<- order(likes.temp$Prior)
  likes.temp$NEVA.Rank<- order(likes.temp$Posterior)
  
  # Need to make predictions from each of these iterations
  mods.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[2]])
  colnames(mods.temp)<- gam.coef
  
  # Get the SDM fits...
  sdm.best<- likes.temp[likes.temp$SDM.Rank<= 900, ]
  sdm.best<- sdm.best[order(sdm.best$SDM.Rank), ]
  sdm.best.params<- mods.temp[sdm.best$Iteration, ] %>%
    summarize_all(., funs(mean, sd), na.rm = T)
  
  # Get the NEVA predictions
  neva.best<- likes.temp[likes.temp$NEVA.Rank<= 900, ]
  neva.best<- neva.best[order(neva.best$NEVA.Rank), ]
  neva.best.params<- mods.temp[neva.best$Iteration, ] %>%
    summarize_all(., funs(mean, sd), na.rm = T)
  
  if(i == 1){
    sdm.result<- data.frame("COMNAME" = spp, "Season" = season, sdm.best.params)
    neva.result<- data.frame("COMNAME" = spp, "Season" = season, neva.best.params)
  } else {
    sdm.temp<- data.frame("COMNAME" = spp, "Season" = season, sdm.best.params)
    sdm.result<- bind_rows(sdm.result, sdm.temp)
    neva.temp<- data.frame("COMNAME" = spp, "Season" = season, neva.best.params)
    neva.result<- bind_rows(neva.result, neva.temp)
  }
}

# Make some plots...
sdm.plot<- sdm.result[,c(1:21)] %>%
  gather(., "Variable", "Mean", -COMNAME, -Season)
sdm.sd<- sdm.result[,c(1, 2, 22:40)] %>%
  gather(., "Variable", "SD", -COMNAME, -Season)
sdm.plot<- bind_cols(sdm.plot, sdm.sd[,-c(1:2)])
sdm.plot$Approach<- rep("SDM", nrow(sdm.plot))

neva.plot<- neva.result[,c(1:21)] %>%
  gather(., "Variable", "Mean", -COMNAME, -Season)
neva.sd<- neva.result[,c(1, 2, 22:40)] %>%
  gather(., "Variable", "SD", -COMNAME, -Season)
neva.plot<- bind_cols(neva.plot, neva.sd[,-c(1:2)])
neva.plot$Approach<- rep("SDM + NEVA", nrow(neva.plot))

params.plot<- bind_rows(sdm.plot, neva.plot)

sdm.plot.fall<- ggplot(na.omit(subset(params.plot, Season == "FALL" & Variable == "X.Intercept._mean")), aes(x = COMNAME, y = Mean, color = Approach)) + 
  scale_color_manual(name = "Approach", values = c("#33a02c", "#1f78b4")) +
  geom_point(position = position_dodge(0.8))+
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width=.2,
                position=position_dodge(0.8)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 

## That didn't really work...what about prediction curves? 
# Results
for(i in seq_along(res.files)){
  
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  
  # Jon's data
  vuln<- paste(jon.df$SENSITIVITY[jon.df$COMNAME == spp], "_", jon.df$SENSITIVITY.CERTAINTY[jon.df$COMNAME == spp], sep = "")
  dir<- paste(jon.df$DIRECTIONAL.EFFECT[jon.df$COMNAME == spp], "_", jon.df$DIRECTIONAL.EFFECT.CERTAINTY[jon.df$COMNAME == spp], sep = "")
  plot.title<- paste("Sens: ", vuln, " and Dir: ", dir, sep = "")
  
  # Model fit -- presence and biomass
  gam.b.temp<- readRDS(paste(out.dir.stem, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(gam.b.temp)$linkinv
  gam.coef<- names(coef(gam.b.temp))
  
  # Ranks of candidate draws...
  # Get ranks for SDM only and for NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  likes.temp$SDM.Rank<- order(likes.temp$Prior)
  likes.temp$NEVA.Rank<- order(likes.temp$Posterior)
  
  # Need to make predictions from each of these iterations
  mods.temp<- data.frame(read_rds(paste(out.dir.stem, res.files[i], sep = ""))[[2]])
  colnames(mods.temp)<- gam.coef
  
  # Get the SDM fits...
  sdm.best<- likes.temp[likes.temp$SDM.Rank<= 900, ]
  sdm.best<- sdm.best[order(sdm.best$SDM.Rank), ]
  sdm.best.params<- mods.temp[sdm.best$Iteration, ] 
  
  # Get the NEVA predictions
  neva.best<- likes.temp[likes.temp$NEVA.Rank<= 900, ]
  neva.best<- neva.best[order(neva.best$NEVA.Rank), ]
  neva.best.params<- mods.temp[neva.best$Iteration, ] 
  
  # Plot the curves...
  # All possible prediction curves, the original smooth, and then the one we have selected...
  # Some prep
  pred.dat.use<- pred.dat[pred.dat$SEASON == season, ]
  rescaled.dat.use<- rescaled.dat[rescaled.dat$Season == season, ]
  
  # Predictor matrix
  pred.mat<- predict(gam.b.temp, newdata = pred.dat.use, type = "lpmatrix")
  
  # All predictions
  pred.all.sdm<- apply(sdm.best.params, 1, function(x) exp(ilink(pred.mat %*% x)))
  pred.all.neva<- apply(neva.best.params, 1, function(x) exp(ilink(pred.mat %*% x)))
  
  # Plotting
  # SST
  want<- 1:500
  sst.all.sdm<- data.frame(pred.all.sdm[want,])
  colnames(sst.all.sdm)<- paste("Mod.", seq(from = 1, to = ncol(sst.all.sdm)), sep = "")
  sst.all.sdm<- sst.all.sdm %>%
    gather(., "Model", "Pred")
  sst.all.sdm$Value<- rep(rescaled.dat.use$SST, length(unique(sst.all.sdm$Model)))
  sst.all.sdm$Parameter<- rep("SST", nrow(sst.all.sdm))
  
  sst.all.neva<- data.frame(pred.all.neva[want,])
  colnames(sst.all.neva)<- paste("Mod.", seq(from = 1, to = ncol(sst.all.neva)), sep = "")
  sst.all.neva<- sst.all.neva %>%
    gather(., "Model", "Pred")
  sst.all.neva$Value<- rep(rescaled.dat.use$SST, length(unique(sst.all.neva$Model)))
  sst.all.neva$Parameter<- rep("SST", nrow(sst.all.neva))
  
  # Mean, Max and Min at each SST value across all models...
  sst.summarized.sdm<- sst.all.sdm %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  sst.summarized.neva<- sst.all.neva %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  sst.dat<- data.frame("Parameter" = rep("SST", 1000), "Value" = rep(rescaled.dat.use$SST, 2), "Pred.Mean" = c(sst.summarized.sdm$mean, sst.summarized.neva$mean), "Pred.Min" = c(sst.summarized.sdm$min, sst.summarized.neva$min), "Pred.Max" = c(sst.summarized.sdm$max, sst.summarized.neva$max), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  sst.out<- ggplot() +
    geom_line(data = sst.dat, aes(x = Value, y = log(Pred.Mean+1), group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("SST") +
    theme_bw() 
  
  sst.out<- sst.out +
    geom_line(data = sst.dat, aes(x = Value, y = log(Pred.Min+1), group = Model, color = Model), lty = "dashed")
  
  sst.out<- sst.out +
    geom_line(data = sst.dat, aes(x = Value, y = log(Pred.Max+1), group = Model, color = Model), lty = "dashed") 
  
  # Depth
  want<- 501:1000
  depth.all.sdm<- data.frame(pred.all.sdm[want,])
  colnames(depth.all.sdm)<- paste("Mod.", seq(from = 1, to = ncol(depth.all.sdm)), sep = "")
  depth.all.sdm<- depth.all.sdm %>%
    gather(., "Model", "Pred")
  depth.all.sdm$Value<- rep(rescaled.dat.use$Depth, length(unique(depth.all.sdm$Model)))
  depth.all.sdm$Parameter<- rep("Depth", nrow(depth.all.sdm))
  
  depth.all.neva<- data.frame(pred.all.neva[want,])
  colnames(depth.all.neva)<- paste("Mod.", seq(from = 1, to = ncol(depth.all.neva)), sep = "")
  depth.all.neva<- depth.all.neva %>%
    gather(., "Model", "Pred")
  depth.all.neva$Value<- rep(rescaled.dat.use$Depth, length(unique(depth.all.neva$Model)))
  depth.all.neva$Parameter<- rep("Depth", nrow(depth.all.neva))
  
  # Mean, Max and Min at each depth value across all models...
  depth.summarized.sdm<- depth.all.sdm %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  depth.summarized.neva<- depth.all.neva %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  depth.dat<- data.frame("Parameter" = rep("Depth", 1000), "Value" = rep(rescaled.dat.use$Depth, 2), "Pred.Mean" = c(depth.summarized.sdm$mean, depth.summarized.neva$mean), "Pred.Min" = c(depth.summarized.sdm$min, depth.summarized.neva$min), "Pred.Max" = c(depth.summarized.sdm$max, depth.summarized.neva$max), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  depth.out<- ggplot() +
    geom_line(data = depth.dat, aes(x = Value, y = log(Pred.Mean+1), group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("Depth") +
    theme_bw() 
  
  depth.out<- depth.out +
    geom_line(data = depth.dat, aes(x = Value, y = log(Pred.Min+1), group = Model, color = Model), lty = "dashed")
  
  depth.out<- depth.out +
    geom_line(data = depth.dat, aes(x = Value, y = log(Pred.Max+1), group = Model, color = Model), lty = "dashed") 
  
  out<- plot_grid(sst.out + theme(legend.position="none"), depth.out + theme(legend.position="none"), nrow = 1, ncol = 2, align = "hv", scale = 1)
  legend<- get_legend(sst.out)
  out<- plot_grid(out, legend, rel_widths = c(3, 0.5))
  title<- ggdraw() + draw_label(plot.title, fontface='bold')
  out<- plot_grid(title, out, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  ggsave(paste(out.dir.stem, tolower(spp), "_", tolower(season), "_PredCurvInt.jpg", sep = ""), out, width = 11, height = 8, dpi = 125, units = "in")
}


######
## Port differences
##### 
fish_avail_func<- function(df) {
  
  if(FALSE){
    row.use<- 1
    data = dat.full$SDM.Diff[[row.use]]
    update = paste(dat.full$COMNAME[[row.use]], dat.full$SEASON[[row.use]], sep = ".")
  }
  
  dat.sp<- data.frame("x" = df$x, "y" = df$y, "z" = df$Projection)
  coordinates(dat.sp)<- ~x+y
  proj4string(dat.sp)<- proj.wgs84
  
  pts.rast<- rasterize(dat.sp, port.foots.stack[[1]], field = "z", fun = mean)
  pts.rast<- raster::resample(pts.rast, port.foots.stack[[1]])
  
  res.mean<- vector(length = raster::nlayers(port.foots.stack), mode = "double")
  res.mean[]<- NA
  res.sd<- vector(length = raster::nlayers(port.foots.stack), mode = "double")
  res.sd[]<- NA
  
  for(i in 1:raster::nlayers(port.foots.stack)) {
    lay.use<- port.foots.stack[[i]]
    
    if(all(is.na(values(lay.use)))){
      res.mean[i]<- NA
      res.sd[i]<- NA
    } else {
      m <- c(0, Inf, 1,  -Inf, 0, 0)
      rclmat <- matrix(m, ncol=3, byrow=TRUE)
      lay.bin<- raster::reclassify(lay.use, rclmat)
      
      # Get coordinates of footprint
      foot.pts <- data.frame(rasterToPoints(lay.bin, function(x) x == 1))
      coordinates(foot.pts)<- ~x+y
      proj4string(foot.pts)<- proj.wgs84
      
      # Okay, we need the projected values at those points and then the proportion of catch at each to use as the weights. 
      # Projected p(presence)
      proj.vals<- raster::extract(pts.rast, foot.pts)
      # Proportion of catch
      proj.weights<- raster::extract(lay.use, foot.pts)
      
      # Mean and SD
      if(all(is.na(proj.vals))) {
        res.mean[i]<- NA
        res.sd[i]<- NA
      } else {
        if(length(proj.vals) == 1){
          res.mean[i]<- mean(proj.vals, na.rm = T)
          res.sd[i]<- NA
        } else {
          res.mean[i]<- mean(proj.vals, na.rm = T)
          res.sd[i]<- sd(proj.vals, na.rm = T)
        }
      }
    }
  }
  
  out<- data.frame("Port" = names(port.foots.stack), "Mean.Avail" = res.mean, "SD.Avail" = res.sd)
  return(out)
}

# Bring in fishing footprints
all.foot.dat<- readRDS("~/GitHub/COCA/Data/VTR fishing footprints by community and gear type 2011-2015.rds")
ports.names<- all.foot.dat$JGS.COMMUNITY

# Spatial projections
proj.wgs84<- "+init=epsg:4326" #WGS84
proj.utm<- "+init=epsg:2960" #UTM 19

# Some work on the names
ports.names<- gsub("\\_(?=[^_]*\\_)", " ", ports.names, perl = TRUE)
ports.names<- gsub(' +', ' ', ports.names)
ports.names<- gsub("_", ".", ports.names)
ports.names<- gsub(" ", "_", ports.names)
ports.names<- gsub("/", "_", ports.names)

# Let's get the lat and long for these areas -- first write it out to do some small edits
if(FALSE){
  write.csv(unique(ports.names), file = "~/Desktop/COCA_FishingPorts_NamesOnly.csv")
  ports.longlat<- read.csv("~/Desktop/COCA_FishingPorts.csv")
  geos.longlats<- geocode(location = paste(ports.longlat$PortName, ports.longlat$PortState), output = "latlon")
  
  # First go....
  ports.longlat$PortLat<- geos.longlats$lat
  ports.longlat$PortLong<- geos.longlats$lon
  
  # Missing after first time
  ports.longlats.nas<- ports.longlat %>%
    filter(., is.na(PortLong))
  ports.longlats.nas2<- geocode(location = paste(ports.longlats.nas$PortName, ports.longlats.nas$PortState), output = "latlon")
  ports.longlats.nas$PortLat<- ports.longlats.nas2$lat
  ports.longlats.nas$PortLong<- ports.longlats.nas2$lon
  
  # Bind....
  ports.longlat<- ports.longlat %>%
    left_join(ports.longlats.nas)
  
  # Another go
  ports.longlats.nas2<- ports.longlat %>%
    filter(., is.na(PortLat))
  ports.longlats.nas3<- geocode(location = paste(ports.longlats.nas2$PortName, ports.longlats.nas2$PortState), output = "latlon")
  ports.longlats.nas2$PortLat<- ports.longlats.nas3$lat
  ports.longlats.nas2$PortLong<- ports.longlats.nas3$lon
  
  # Bind...
  ports.longlat<- ports.longlat %>%
    left_join(ports.longlats.nas2)
  
  # Still some missing...write it out and add em manually
  write.csv(ports.longlat, file = "~/Desktop/COCA_FishingPorts.csv")
  ports.longlat<- read.csv("~/Desktop/COCA_FishingPorts.csv")
}

ports.longlat<- read.csv("~/Desktop/NormalVoting_BiomassIncPresNoExposure_09222018/COCA_FishingPorts.csv")

# Fishing port footprints -- gear type specific
gear.types<- all.foot.dat$COST_ID
gear.types<- ifelse(gear.types == 1, "Dredge",
                    ifelse(gear.types == 2, "Gillnet",
                           ifelse(gear.types == 3, "Longline",
                                  ifelse(gear.types == 4, "Pot/Trap",
                                         ifelse(gear.types == 5, "Purse/Seine",
                                                ifelse(gear.types == 6, "Trawl", "Other"))))))
port.foot.names<- paste(ports.names, gear.types, sep = "-")
ports.all.foots<- all.foot.dat$JGS.NOAA.SAFE.COMMUNITY.GEAR.FOOTPRINTS
names(ports.all.foots)<- port.foot.names

# Get proportion layer we want for each port-gear type
port.foots<- unlist(lapply(ports.all.foots, "[", 3))
port.foots.stack<- raster::stack(port.foots)
names1<- names(port.foots.stack)
rast.ind<- nlayers(port.foots.stack)

# Also need an "all gear" option and a max distance option
ports.only<- str_replace_all(names(port.foots.stack), c(".Pot.Trap.JGS.SAFE.PROPORTION" = "", ".Other.JGS.SAFE.PROPORTION" = "", ".Gillnet.JGS.SAFE.PROPORTION" = "", ".Trawl.JGS.SAFE.PROPORTION" = "", ".Dredge.JGS.SAFE.PROPORTION" = "", ".Purse.Seine.JGS.SAFE.PROPORTION" = "", ".Longline.JGS.SAFE.PROPORTION" = ""))
ports.unique<- unique(ports.only)

all.stack<- stack()

for(i in 1:length(ports.unique)){
  port.use<- ports.unique[i]
  port.foots.ind<- which(grepl(port.use, names(port.foots.stack)), arr.ind = T)
  stack.use<- port.foots.stack[[port.foots.ind]]
  if(nlayers(stack.use) == 1 | all(is.na(getValues(stack.use)))){
    all.gear.out<- stack.use[[1]]
  } else {
    all.gear.out<- sum(stack.use, na.rm = T)
  }
  all.stack<- raster::stack(all.stack, all.gear.out)
}

# Combine
names.all<- c(names1, paste(ports.unique, ".All.JGS.SAFE.PROPORTION", sep = ""))
check.stack<- raster::stack(port.foots.stack, all.stack)
names(check.stack)<- names.all
port.foots.stack<- check.stack

# Now max distance option
dist.fac<- 1.5

# Similar loop, but we don't want to do this for the "All gear" scenario...
port.foots.stack.noall<- port.foots.stack[[-which(grepl(".All.JGS.SAFE.PROPORTION", names(port.foots.stack)))]]

# Empty stack for results
stack.maxd<- stack()

# Loop
for(i in 1:nlayers(port.foots.stack.noall)){
  lay.use<- port.foots.stack.noall[[i]]
  
  if(all(is.na(getValues(lay.use)))){
    max.d<- lay.use
    max.d.fac<- max.d
    stack.maxd<- raster::stack(stack.maxd, max.d, max.d.fac)
    next()
  } else {
    port.name.use<- str_replace_all(names(lay.use), c(".Pot.Trap.JGS.SAFE.PROPORTION" = "", ".Other.JGS.SAFE.PROPORTION" = "", ".Gillnet.JGS.SAFE.PROPORTION" = "", ".Trawl.JGS.SAFE.PROPORTION" = "", ".Dredge.JGS.SAFE.PROPORTION" = "", ".Purse.Seine.JGS.SAFE.PROPORTION" = "", ".Longline.JGS.SAFE.PROPORTION" = ""))
    pt.use<- ports.longlat[ports.longlat$JGS.PortName == port.name.use, ]
    coordinates(pt.use)<- ~PortLong+PortLat
    proj4string(pt.use)<- proj.wgs84
    
    # Get distance from port, then get maximum distance given fishing for max dist and maxdist * 1.5
    dist.rast<- distanceFromPoints(lay.use, pt.use)
    max.d<- dist.rast
    max.d.maxval<- maxValue(mask(dist.rast, lay.use))
    max.d[max.d >= max.d.maxval]<- NA
    
    max.d.fac<- dist.rast
    max.d.fac[max.d.fac >= max.d.maxval*dist.fac]<- NA
    
    stack.maxd<- raster::stack(stack.maxd, max.d, max.d.fac)
  }
}

# Combine
names.stackmaxd<- paste(rep(names(port.foots.stack.noall), each = 2), c("MaxD.JGS.SAFE.PROPORTION", "1.5xMaxD.JGS.SAFE.PROPORTION"), sep = "")
names(stack.maxd)<- names.stackmaxd
names.all<- c(names(port.foots.stack), names.stackmaxd)
check.stack<- raster::stack(port.foots.stack, stack.maxd)
names(check.stack)<- names.all
port.foots.stack<- check.stack

# Read in results and apply fish availability function to each of the projections datasets
results<- read_rds(paste(out.dir.stem, "atlantic surfclam", "SDMPredictions.rds", sep = "")) %>%
  mutate(., "Fish.Availability" = purrr::map(Projections, fish_avail_func))

saveRDS(results, paste(out.dir.stem, "atlantic surfclam", "FishAvailability.rds", sep = ""))

# Ideally, we want species-season-port and then all the changes as columns...
results<- readRDS(paste(out.dir.stem, "atlantic surfclam", "FishAvailability.rds", sep = ""))
results$spp.season<- paste(results$COMNAME, results$SEASON, sep = ".")
spp.season<- paste(rep(unique(results$COMNAME), each = 1), rep("SUMMER"), sep = ".")

for(i in seq_along(spp.season)){
  dat.temp<- results[results$spp.season == spp.season[i],]
  
  if(nrow(dat.temp) == 0){
    print(paste(spp.season[i], " is done!", sep = " "))
    next()
  }
  
  dat.temp<- dat.temp %>% 
    dplyr::select(., COMNAME, SEASON, Proj.Class, Fish.Availability) %>%
    unnest(Fish.Availability) %>%
    gather(., "Stat", "Value", -COMNAME, -SEASON, -Proj.Class, -Port) %>%
    unite(Port.Stat, Port, Stat) %>%
    spread(Proj.Class, Value) %>%
    separate(Port.Stat, c("Port", "Stat"), sep = "_[^_]*$")
  dat.temp$Stat<- rep(c("Mean", "SD"), length.out = nrow(dat.temp))
  
  if(i == 1){
    result<- dat.temp
    print(paste(spp.season[i], " is done!", sep = " "))
  } else {
    result<- bind_rows(result, dat.temp)
    print(paste(spp.season[i], "is done!", sep = " "))
  }
}

# Group by species and Port, take column means
out<- result %>% 
  group_by(COMNAME, Port, Stat) %>%
  summarize_if(is.double, mean, na.rm = T) %>%
  arrange(., COMNAME, Stat, Port)

# Separate Port and Gear Type
port.data.out<- out
names(port.data.out)[2]<- "Port_GearType"

# First Do gear types...
port.data.out$Gear<- ifelse(grepl(".Pot.Trap", port.data.out$Port_GearType), "Pot.Trap",
                    ifelse(grepl(".Gillnet", port.data.out$Port_GearType), "Gillnet",
                           ifelse(grepl(".Longline", port.data.out$Port_GearType), "Longline",
                                  ifelse(grepl(".Dredge", port.data.out$Port_GearType), "Dredge",
                                         ifelse(grepl(".Purse.Seine", port.data.out$Port_GearType), "Purse.Seine",
                                                ifelse(grepl(".Trawl", port.data.out$Port_GearType), "Trawl",
                                                       ifelse(grepl(".All", port.data.out$Port_GearType), "All", "Other")))))))
port.data.out$Footprint<- ifelse(str_detect(port.data.out$Port_GearType, "1.5xMaxD.JGS.SAFE.PROPORTION"), "1.5xMaxD", "Regular")
port.data.out$Footprint<- ifelse(str_detect(port.data.out$Port_GearType, "MaxD.JGS.SAFE.PROPORTION") & port.data.out$Footprint == "Regular", "MaxD", port.data.out$Footprint)

port.data.out$Port.OnlyTemp<- str_replace_all(port.data.out$Port_GearType, c(".JGS.SAFE.PROPORTION"), "") 
port.data.out$Port.OnlyTemp2<- str_replace_all(port.data.out$Port.OnlyTemp, c("1.5xMaxD"), "")
port.data.out$Port.OnlyTemp3<- str_replace_all(port.data.out$Port.OnlyTemp2, c("MaxD"), "")
port.data.out$Port.Only<- str_replace_all(port.data.out$Port.OnlyTemp3, c(".Pot.Trap" = "", ".Dredge" = "", ".Longline" = "", ".All" = "", ".Other" = "", ".Gillnet" = "", ".Trawl" = "", ".Purse.Seine" = ""))

# Community and State
port.data.out$Community.Only<- unlist(lapply(strsplit(sub("\\.", "*", port.data.out$Port.Only), "\\*"), "[", 1))
port.data.out$State.Only<- unlist(lapply(strsplit(sub("\\.", "*", port.data.out$Port.Only), "\\*"), "[", 2))

port.data.out<- port.data.out %>%
  ungroup() %>%
  dplyr::select(., COMNAME, Port_GearType, Community.Only, State.Only, Gear, Footprint, colnames(port.data.out)[3:24])

# Merge with Brad's port names
merge.table<- read_csv("~/GitHub/COCA/Data/JGS_vs_CFDERRS_PortNamesMatchTable.csv")
merge.table$Merge.Col<- gsub(" ", "_", merge.table$`JGS FULL NAME`)

port.data.out$Merge.Col<- paste(gsub("[.]", "_", gsub(" ", "_", port.data.out$Community.Only)), port.data.out$State.Only, sep = "_")
port.data.out$Econ.Port.Name<- merge.table$`BRAD PORT_NAME_STATE`[match(port.data.out$Merge.Col, merge.table$Merge.Col)]
port.data.out<- port.data.out %>%
  dplyr::select(., COMNAME, Port_GearType, Community.Only, State.Only, Econ.Port.Name, Gear, Footprint, colnames(port.data.out)[7:28])

names(port.data.out)<- gsub(x = names(port.data.out), pattern = "\\.", replacement = "_")
port.data.out$Gear<- gsub("[.]", "_", port.data.out$Gear)
write.csv(port.data.out, "~/Desktop/AtlanticSurfclamPortData09222018.csv")

clams<- read_csv("~/Desktop/AtlanticSurfclamPortData09222018.csv")
all<- read_csv("~/Desktop/PortData09122018.csv") %>%
  bind_rows(., clams) %>%
  arrange(., COMNAME, Econ_Port_Name, Footprint, Stat)
write.csv(all, "~/Desktop/AllPortData09222018.csv")

brad<- all %>%
  filter(., Stat == "Mean")
write.csv(brad, "~/Desktop/BradPortData09222018.csv")

km<- brad %>%
  filter(., grepl("All", Gear)) %>%
  dplyr::select(., COMNAME, Community_Only, State_Only, Gear, Footprint, Baseline_combo_b, Baseline_sdm_b, Baseline_sdm_p, Future_mean_combo_b, Future_mean_sdm_b, Future_mean_sdm_p, Future_mean_diff_combo_b, Future_mean_diff_sdm_b, Future_mean_diff_sdm_p)
write.csv(km, "~/Desktop/KMPortData09222018.csv")

#####
## Shelfwide and regional changes
#####
results<- read_rds(paste(out.dir.stem, "SDMPredictions.rds", sep = "")) # This should have everything we need. Let's get the model fit, too....
if(FALSE){
  dat.full<- dat.full.keep
}
mod.results<- read.csv(paste(out.dir.stem, "mod.results.csv", sep = ""))
dat.full<- results %>%
  left_join(., mod.results, by = c("COMNAME", "SEASON")) %>%
  dplyr::select(., -X)
dat.full$AUC.SDM[is.na(dat.full$AUC.SDM)]<- 0

dat.sub<- dat.full %>%
  filter(., AUC.SDM >= 0.70)

dat.full<- dat.sub

# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

# NELME
nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
st_crs(nelme)<- "+init=epsg:4326"
nelme.sp<- as(nelme, "Spatial")

# GoM
gom<- st_read("~/Dropbox/Andrew/Work/GMRI/AllGIS/PhysioRegions_Maine/PhysioRegions_wgs84.shp")
st_crs(gom)<- "+init=epsg:4326"
gom<- gom[!gom$Region == "Seamount",] %>%
  st_union() %>%
  st_sf()
gom.sp<- as(st_zm(gom), "Spatial")
gom.sp<- spTransform(gom.sp, proj.utm)
gom.sp.wgs<-  spTransform(gom.sp, proj4string(nelme.sp))

# Buffer it a bit
gom.buff<- gBuffer(gom.sp, width = 7500)
gom.buff<- spTransform(gom.buff, proj4string(nelme.sp))

# Southern regions
south<- erase(nelme.sp, gom.buff)

# Still a bit remaining...custom box to get rid of the rest of it
# Coordinates
ow<- data.frame("x" = c(-71, -71, -67, -67), "y" = c(42, 46, 46, 42))

# Convert coordinates to Spatial Polygons
ow.p<- Polygon(ow)
ow.ps<- Polygons(list(ow.p), 1)
ow.sp<- SpatialPolygons(list(ow.ps))
proj4string(ow.sp)<- proj4string(nelme.sp)
south2<- erase(south, ow.sp)
proj4string(south2)<- proj4string(nelme.sp)

# Overlay func
overlay_func<- function(df, region, proj.use = proj4string(nelme.sp)){
  dat.use<- data.frame(df)
  pts.temp<- dat.use
  coordinates(pts.temp)<- ~x+y
  proj4string(pts.temp)<- proj.use
  
  switch(region,
         NELME = mean(dat.use[,3], na.rm = T),
         GOM = mean(data.frame(pts.temp[!is.na(over(pts.temp, as(gom.sp.wgs, "SpatialPolygons"))),])[,3], na.rm = T),
         South = mean(data.frame(pts.temp[!is.na(over(pts.temp, as(south2, "SpatialPolygons"))),])[,3], na.rm = T))
}

preds.df.sub<- dat.full %>%
  mutate(., "NELME.Mean.Change" = purrr::map2(Projections, list("NELME"), possibly(overlay_func, NA)),
         "GOM.Mean.Change" = purrr::map2(Projections, list("GOM"), possibly(overlay_func, NA)),
         "South.Mean.Change" = purrr::map2(Projections, list("South"), possibly(overlay_func, NA)))


## Alright, we are now after a species - scenario - season - region - mean change dataframe...
res<- preds.df.sub %>%
  dplyr::select(., COMNAME, SEASON, Proj.Class, NELME.Mean.Change, GOM.Mean.Change, South.Mean.Change) %>%
  gather(., "Region", "Change", -COMNAME, -SEASON, -Proj.Class) 
res$Change<- as.numeric(unlist(res$Change))

res.p.vec<- c("Future_mean_diff.sdm.p", "Future_cold_diff.sdm.p", "Future_warm_diff.sdm.p")
res.b.vec<- c("Future_mean_diff.combo.b", "Future_cold_diff.combo.b", "Future_warm_diff.combo.b")

res.plot<- res %>%
  dplyr::filter(., Proj.Class %in% c(res.p.vec, res.b.vec))
res.plot$Region<- gsub(".Mean.Change", "", res.plot$Region)

# Lets add a functional group column...
# Merge with functional groups....
func.groups<- read.csv("~/Dropbox/Andrew/Work/GMRI/COCA/Data/JHareSppFunctionalGroup.csv")
func.groups$COMNAME<- toupper(func.groups$COMNAME)

res.plot<- res.plot %>%
  left_join(., func.groups, by = "COMNAME")
res.plot<- res.plot[!is.na(res.plot$Functional.Group),]
res.plot$Functional.Group<- factor(res.plot$Functional.Group, levels = c("Groundfish", "Pelagic", "Coastal", "Invertebrates", "Diadromous", "Elasmobranch"))
res.plot<- res.plot %>%
  dplyr::arrange(., COMNAME, Functional.Group, SEASON, Proj.Class, Region, Change)
res.plot$COMNAME<- factor(res.plot$COMNAME, levels = unique(res.plot$COMNAME))
res.plot$Model<- ifelse(grepl("combo", res.plot$Proj.Class), "Biomass", "Presence")
res.plot$SEASON<- factor(res.plot$SEASON, levels = c("FALL", "SPRING"))
res.plot$Proj.Class<- factor(res.plot$Proj.Class, levels = c(res.p.vec, res.b.vec))
res.plot$Region<- factor(res.plot$Region, levels = c("NELME", "GOM", "South"))
res.plot$Model<- factor(res.plot$Model, levels = c("Presence", "Biomass"))

res.plot.nelme<- res.plot %>%
  dplyr::filter(., Region == "NELME")
res.plot.gom<- res.plot %>%
  dplyr::filter(., Region == "GOM")
res.plot.south<- res.plot %>%
  dplyr::filter(., Region == "South")

res.plot.all<- list(res.plot.nelme, res.plot.gom, res.plot.south)
names(res.plot.all)<- c("NELME", "GoM", "South")

for(i in seq_along(res.plot.all)){
  df<- data.frame(res.plot.all[[i]])
  df.null<- cbind(expand.grid(COMNAME = levels(res.plot.nelme$COMNAME), SEASON = unique(res.plot.nelme$SEASON), Proj.Class = unique(res.plot.nelme$Proj.Class), Model = levels(res.plot.nelme$Model), Region = levels(res.plot.nelme$Region), Change = NA))
  df.null<- df.null %>%
    left_join(., func.groups, by = "COMNAME")
  df<- rbind(df[,], df.null)
  df$duplicated<- paste(df$COMNAME, df$SEASON, df$Proj.Class, df$Functional.Group)
  df<- df[!duplicated(df$duplicated),] %>%
    arrange(., COMNAME, SEASON)
  
  for(j in seq_along(levels(df$SEASON))){
    dat.use<- df %>%
      dplyr::filter(., SEASON == levels(df$SEASON)[j])
    
    dodge <- position_dodge(width = 1)
    
    dat.use.df<- dat.use %>%
      tidyr::complete(COMNAME, SEASON)
    dat.use.df$COMNAME<- str_to_title(dat.use.df$COMNAME)
    dat.use.df$COMNAME<- factor(dat.use.df$COMNAME, levels = rev(unique(dat.use.df$COMNAME)))
    
    dat.use.df<- dat.use.df %>%
      drop_na(Change)
    dat.use.df$Count<- ifelse(dat.use.df$Change == 0, 0, ifelse(dat.use.df$Change > 0, 2, -1))
    grouped.df<- dat.use.df %>%
      group_by(., COMNAME) %>%
      dplyr::summarize(., "Plot.Group" = sum(Count))
    
    # Join
    dat.use.df<- dat.use.df %>%
      left_join(., grouped.df, by = "COMNAME")
    dat.use.df$Plot.Group<- factor(dat.use.df$Plot.Group, levels = rev(c(-2, -1, 1, 0, 2, 4)))
    dat.use.df<- dat.use.df %>%
      arrange(., Plot.Group, COMNAME, SEASON)
    dat.use.df$COMNAME.Plot<- factor(dat.use.df$COMNAME, levels = unique(dat.use.df$COMNAME))
    
    means.df<- dat.use.df %>%
      dplyr::filter(., grepl("mean", Proj.Class))
    scenarios.df<- dat.use.df %>%
      dplyr::filter(., !grepl("mean", Proj.Class))
    scenarios.df$Scenario<- ifelse(grepl("cold", scenarios.df$Proj.Class), "Cold", "Warm")
    scenarios.df$Scenario<- factor(scenarios.df$Scenario, levels = c("Cold", "Warm"))
    
    plot.means<- ggplot(data = means.df, aes(COMNAME.Plot, Change, fill = Model)) + 
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.6)) +
      scale_fill_manual("Response", values = c("#cccccc", "#969696")) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      theme(text = element_text(size = 12)) +
      coord_flip() +
      facet_wrap(~Functional.Group, scales = "free") +
      ggtitle(paste(names(res.plot.all)[i], levels(df$SEASON)[j], sep = " "))
    plot.error<- plot.means +
      geom_point(data = scenarios.df, aes(COMNAME.Plot, Change, color = Scenario, group = Model), position = position_dodge(width = 0.6)) +
      scale_color_manual("Climate Scenario", values = c("#3182bd", "#de2d26")) +
      coord_flip() +
      facet_wrap(~Functional.Group, scales = "free")
    
    ggplot2::ggsave(filename = paste(out.dir.stem, names(res.plot.all)[i], levels(df$SEASON)[j], ".jpg", sep = ""), plot = plot.error, width = 11, height = 8, units = "in")
  }
  
  print(paste(names(res.plot.all)[i], levels(df$SEASON)[j], "is done!", sep = " "))
}

#####
## Change weighted by landings
#####
port.diffs<- read.csv("~/Desktop/PortData09042018.csv")
port.dat.l<- port.diffs %>%
  dplyr::filter(grepl("All", Gear)) %>%
  group_by(., COMNAME, Port_Long) %>%
  summarize_if(is.double, mean, na.rm = T) %>% 
  drop_na() 
port.dat.l$Port.long<- gsub(".All.JGS.SAFE.PROPORTION", "", port.dat.l$Port_Long)

# Get locations
ports.only<- data.frame(port.dat.l[!duplicated(port.dat.l$Port.long),])
ports.only$Address<- as.character(gsub("[_]", " ", gsub("[.]", " ", tolower(ports.only$Port.long))))
port.location<- mutate_geocode(ports.only, Address)
port.location<- port.location %>%
  dplyr::select(., Port.long, lon, lat)

port.dat.l<- port.dat.l %>%
  left_join(., port.location, by = "Port.long")

# Need to add importance and x/y...
landings.dat<- read.csv("~/GitHub/COCA/Data/landport1115.csv")

# From this wide datafile, we want to add the value where species name matches column names of landings.dat and port name matches port row.
column.indices<- match(gsub(" ", ".", port.dat.l$COMNAME), colnames(landings.dat))
row.indices<- match(port.dat.l$Port.long, landings.dat$Port.long)
landings.dat$filler<- rep(NA, nrow(landings.dat))

# Adjust NA columns (no match in species name)
column.indices[is.na(column.indices)]<- ncol(landings.dat)

# Add the landings and calculate proportion of landings...by port?
port.dat.l$landings<- as.numeric(landings.dat[cbind(row.indices, column.indices)])

port.dat.weighted.diff<- port.dat.l %>%
  dplyr::group_by(., Port.long) %>%
  dplyr::mutate(weighted.importance = round(landings/sum(landings, na.rm = T), 3),
                weighted.difference.p = Future_mean_diff_sdm_p*weighted.importance,
                weighted.difference.b = Future_mean_diff_combo_b*weighted.importance) %>%
  dplyr::select(., COMNAME, Port_Long, lon, lat, landings, weighted.importance, weighted.difference.p, weighted.difference.b)

# Write this out
write.csv(data.frame(port.dat.weighted.diff), file = paste(out.dir.stem, "SpeciesPortGearTypeChanges.csv"))

# Calculate port specific sum
port.aggregated.weighted.difference<- port.dat.weighted.diff %>%
  dplyr::group_by(Port.long, lon, lat) %>%
  dplyr::mutate(., TotalChange.P = sum(weighted.difference.p, na.rm = T),
                TotalChange.B = sum(weighted.difference.b, na.rm = T)) %>%
  data.frame()

# Write this out and return it
write.csv(data.frame(port.aggregated.weighted.difference), file = paste(out.dir.stem, "PortAggregatedImportanceWeightedChange.csv"))

# Plots for focal ports
focal.ports<- c("STONINGTON.ME", "PORTLAND.ME", "NEW_BEDFORD.MA", "POINT_JUDITH.RI")
port.imp.plot<- port.aggregated.weighted.difference %>%
  filter(., Port.long %in% focal.ports)
port.imp.plot<- port.imp.plot %>%
  dplyr::select(COMNAME, Port_Long, TotalChange.P, TotalChange.B) %>%
  gather(., "Response", "TotalChange", -COMNAME, -Port_Long)

port.imp.plot.out<- ggplot(data = port.imp.plot, aes(x = Port_Long, y = TotalChange, fill = Response, group = Response)) +
  geom_bar(stat = "identity", position = "dodge")  +
  ylab("Proportional Landings\n Weighted Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size = 18))
port.imp.plot.out.b<- ggplot(data = port.imp.plot, aes(x = Port.Name, y = TotalChange.B)) +
  geom_bar(stat = "identity", position = "dodge")  +
  ylab("Proportional Landings\n Weighted Change in P(presence)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size = 18))

png(file = paste(out.dir.stem, "FocalPorts_AverageChangeWeightedByImp.png", sep = ""), width = 12, height = 8, units = "in", res = 250)
plot(port.imp.plot.out)
dev.off()




#####
## Matching up with Jon's designations
#####



