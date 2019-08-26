### One species and One port footprint examples
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

temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

#####
## Get some data for one species
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")
spp.example<- c("LONGFIN SQUID", "AMERICAN LOBSTER")

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

#####
## Fit the model
dat.use<- dat.train[2,]
season.use<- as.character(dat.use$SEASON[[1]])
gam.mod0<- gam(PRESENCE.ABUNDANCE ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)

#####
## Predictions
# SDM
season<- "SPRING"
base.map<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = predict.gam(gam.mod0, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "response"))
fut.map<- data.frame("x" = fut.preds$Data[[match(season, fut.preds$SEASON)]]$x, "y" = fut.preds$Data[[match(season, fut.preds$SEASON)]]$y, "pred" = predict.gam(gam.mod0, newdata = fut.preds$Data[[match(season, fut.preds$SEASON)]], type = "response"))

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

# Baseline -- SDM first and then SDM and NEVA
data.use<- base.map
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

# Discrete scale
pred.df.use$breaks<- cut(pred.df.use$z, 
                         breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                         labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
pred.df.use<- pred.df.use %>%
  drop_na(., breaks)

plot.base.sdm<- ggplot() + 
  geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
  scale_fill_viridis(option = "viridis", name = "Probability of\n presence", na.value = "white", limits = c(0,1)) +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10))

# Add star for Port location
plot.base.sdm<- plot.base.sdm +
  geom_point(aes(x = -70.7626, y = 43.0718), shape = 21, fill = "red", color = "black", size = 3)

# Add footprint
focal.ports<- "PORTSMOUTH_NH"
all.foot.dat<- readRDS("./Data/VTR fishing footprints by community and gear type 2011-2015.rds")
all.foot.dat<- all.foot.dat[all.foot.dat$JGS.COMMUNITY == focal.ports,]
ports.names<- all.foot.dat$JGS.COMMUNITY

gear.types<- all.foot.dat$COST_ID
gear.types<- ifelse(gear.types == 1, "Dredge",
                    ifelse(gear.types == 2, "Gillnet",
                           ifelse(gear.types == 3, "Longline",
                                  ifelse(gear.types == 4, "Pot/Trap",
                                         ifelse(gear.types == 5, "Purse/Seine",
                                                ifelse(gear.types == 6, "Trawl", "Other"))))))
port.foot.names<- paste(ports.names, gear.types, sep = ".")
ports.all.foots<-all.foot.dat$JGS.NOAA.SAFE.COMMUNITY.GEAR.FOOTPRINTS
names(ports.all.foots)<- port.foot.names

# Create a layer of all gear types for each port, add it to our list
unique.ports<- unique(ports.names)

# Add all gear type by port layer to the existing stack
temp.list<- unlist(lapply(ports.all.foots, "[", 3))
start<- length(temp.list)

if(!is.null(focal.ports)) {
  unique.ports<- focal.ports
  
  # Add all gear type by port layer to the existing stack
  temp.list<- unlist(lapply(ports.all.foots, "[", 3))
  start<- length(temp.list)
} else {
  unique.ports<- unique(ports.names)
  
  # Add all gear type by port layer to the existing stack
  temp.list<- unlist(lapply(ports.all.foots, "[", 3))
  start<- length(temp.list)
}

for(i in 1:length(unique.ports)) {
  port.use<- unique.ports[i]
  
  port.dat.all<- ports.all.foots[grepl(port.use, names(ports.all.foots))]
  all.layers<- raster::stack(unlist(lapply(port.dat.all, "[", 3)))
  all.gear.layer<- sum(raster::stack(all.layers), na.rm = T)
  temp.list[[start+i]]<- all.gear.layer
  names(temp.list)[[start+i]]<- paste(port.use, "All.Gear", sep = ".")
  print(paste(port.use, " is done!", sep = ""))
}

# Out
port.foots<- temp.list

# Getting all gear
port.rast.dat<- raster::as.data.frame(port.foots[[7]], xy = TRUE)
names(port.rast.dat)[3]<- "z"

plot.out<- ggplot() + 
  geom_raster(data = port.rast.dat, aes(x = x, y = y, fill = z), show.legend = FALSE) +
  scale_fill_viridis(option = "magma", name = "Proportion of Kept Catch", na.value = "gray80") +
  geom_map(data = us.states.f, map = us.states.f,
           aes(x = long, y = lat, map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(x = long, y = lat, map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  ggtitle(paste(port.use)) +
  theme(panel.background = element_rect(fill = "gray80", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black")) 

# Want this extent
# First, make all values the same. Either do
focal.rast<- port.foots[[7]]

# Next, convert to polygons (you need to have package 'rgeos' installed for this to work)
focal.poly<- rasterToPolygons(focal.rast, dissolve=TRUE)

# Fortify
focal.poly.df<- fortify(focal.poly)

plot.sdm.out<- plot.base.sdm +
  geom_map(data = focal.poly.df, map = focal.poly.df, aes(map_id = id, group = group), fill = NA, color = "red", size = 1)

