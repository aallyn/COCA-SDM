######
## Footprint weirdness
######
library(tidyverse)
library(raster)
library(ggmap)

# Noticed something strange when bringing together the model projections with the landings data...
sdm.path<- "/Volumes/Shared/Research/COCA-conf/SDM and CFDERS Integration/Data/SDM Projections/"
landings.file<- "/Volumes/Shared/Research/COCA-conf/SDM and CFDERS Integration/Processed Summaries/FocalComm.Spp.Gearsummary.csv"
focal.comms<- c("STONINGTON_ME", "PORTLAND_ME", "NEW BEDFORD_MA", "POINT JUDITH_RI")
out.path<- "/Volumes/Shared/Research/COCA-conf/SDM and CFDERS Integration/Processed Summaries/"

# For any of these, same process for dealing with model results
mod.stats<- read_csv(paste(sdm.path, "mod.results.csv", sep = "")) %>%
  dplyr::select(-X1) %>%
  group_by(COMNAME) %>%
  summarize_if(., is.numeric, mean, na.rm = TRUE)
colnames(mod.stats)[1]<- "CommonName"

# First, determine what file is being used for weights CFDERS "FocalComm.Spp.Gear" summary file or "Comm.Spp" summary file or GAR summary file
file.use<- basename(landings.file)

# Read in species changes and remove "All" gear 
comm.diffs<- read_csv(paste(sdm.path, "EcoToEconPortData03032019.csv", sep = "")) %>%
  dplyr::select(., -X1) %>%
  filter(., Community %in% focal.comms & Gear != "All") %>%
  drop_na(CFDERSCommonName, Value)
colnames(comm.diffs)[10]<- "ProjectionValue"

# Merge over model results
comm.diffs<- comm.diffs %>%
  left_join(., mod.stats)

# For each species - community - gear - projection scenario there should be three footprints...
foot.check<- comm.diffs %>%
  group_by(., CommonName, Community, Gear, ProjectionScenario) %>%
  summarize(.,
            "DistinctFoots" = n_distinct(Footprint))
summary(foot.check)

# Clearly not the case...
foots.miss<- foot.check %>%
  filter(., DistinctFoots < 3)

unique(foots.miss$CommonName) # Nothing jumps out here...
unique(foots.miss$ProjectionScenario) # Only happening in percent differences...There's a chance this got fixed with the new code...but just to see...

# What do these footprints look like
# Bring in fishing footprints
all.foot.dat<- readRDS("~/GitHub/COCA/Data/VTR fishing footprints by community and gear type 2011-2015.rds")
ports.names<- all.foot.dat$JGS.COMMUNITY

# Spatial projections
proj.wgs84<- "+init=epsg:4326" #WGS84
proj.utm<- "+init=epsg:2960" #UTM 19

# Some work on the names
ports.names<- gsub("\\_(?=[^_]*\\_)", " ", ports.names, perl = TRUE)
ports.names<- gsub(' +', ' ', ports.names)
ports.names<- gsub("/", " ", ports.names)

port.and.state<- strsplit(ports.names, split = "_")

ports.geo<- data.frame("PortName" = unlist(lapply(port.and.state, "[", 1)), "PortState" = unlist(lapply(port.and.state, "[", 2)))

# Let's get the lat and long for these areas -- first write it out to do some small edits
google.api<- "AIzaSyDqGwT79r3odaMl0Hksx9GZhHmqe37KSEQ"
register_google(key = google.api)
geos.longlats<- geocode(location = unique(paste(ports.geo$PortName, ports.geo$PortState)), output = "latlon")
geo.merge<- data.frame("PortMerge" = unique(paste(ports.geo$PortName, ports.geo$PortState)), "Long" = geos.longlats$lon, "Lat" = geos.longlats$lat)

ports.geo<- ports.geo %>%
  mutate(., PortMerge = paste(PortName, PortState)) %>%
  left_join(., geo.merge)

ports.geo$JGSCommunity<- all.foot.dat$JGS.COMMUNITY

# Fishing port footprints -- gear type specific
gear.types<- all.foot.dat$COST_ID
gear.types<- ifelse(gear.types == 1, "Dredge",
                    ifelse(gear.types == 2, "Gillnet",
                           ifelse(gear.types == 3, "Longline",
                                  ifelse(gear.types == 4, "Pot/Trap",
                                         ifelse(gear.types == 5, "Purse/Seine",
                                                ifelse(gear.types == 6, "Trawl", "Other"))))))
port.foot.names<- paste(ports.names, gear.types, sep = "-")

# Safe vs. unsafe
unsafe<- TRUE
if(unsafe){
  ports.all.foots<- all.foot.dat$JGS.COMMUNITY.GEAR.FOOTPRINTS
} else {
  ports.all.foots<- all.foot.dat$JGS.NOAA.SAFE.COMMUNITY.GEAR.FOOTPRINTS
  
}
ports.all.foots<- all.foot.dat$JGS.COMMUNITY.GEAR.FOOTPRINTS
names(ports.all.foots)<- port.foot.names

# Get proportion layer we want for each port-gear type
port.foots<- unlist(lapply(ports.all.foots, "[", 3))
port.foots.stack<- raster::stack(port.foots) # 476 layers
names1<- names(port.foots.stack)
rast.ind<- nlayers(port.foots.stack)

# Also need an "all gear" option and a max distance option
if(unsafe){
  ports.only<- str_replace_all(names(port.foots.stack), c(".Pot.Trap.JGS.PROPORTION" = "", ".Other.JGS.PROPORTION" = "", ".Gillnet.JGS.PROPORTION" = "", ".Trawl.JGS.PROPORTION" = "", ".Dredge.JGS.PROPORTION" = "", ".Purse.Seine.JGS.PROPORTION" = "", ".Longline.JGS.PROPORTION" = ""))
} else {
  ports.only<- str_replace_all(names(port.foots.stack), c(".Pot.Trap.JGS.SAFE.PROPORTION" = "", ".Other.JGS.SAFE.PROPORTION" = "", ".Gillnet.JGS.SAFE.PROPORTION" = "", ".Trawl.JGS.SAFE.PROPORTION" = "", ".Dredge.JGS.SAFE.PROPORTION" = "", ".Purse.Seine.JGS.SAFE.PROPORTION" = "", ".Longline.JGS.SAFE.PROPORTION" = ""))
}

ports.unique<- unique(ports.only) # 126 ports

all.stack<- stack()

for(i in 1:length(ports.unique)){
  port.use<- ports.unique[i]
  port.foots.ind<- which(grepl(port.use, names(port.foots.stack)), arr.ind = T)
  stack.use<- port.foots.stack[[port.foots.ind]]
  if(nlayers(stack.use) == 1){
    all.gear.out<- stack.use[[1]]
  } else {
    all.gear.out<- calc(stack.use, sum, na.rm = T)
  }
  all.stack<- raster::stack(all.stack, all.gear.out)
}

# Combine
if(unsafe){
  names.all<- c(names1, paste(ports.unique, ".All.JGS.PROPORTION", sep = ""))
} else {
  names.all<- c(names1, paste(ports.unique, ".All.JGS.SAFE.PROPORTION", sep = ""))
}
check.stack<- raster::stack(port.foots.stack, all.stack)
names(check.stack)<- names.all
port.foots.stack<- check.stack

# Now max distance option
dist.fac<- 1.5

# Similar loop, but we don't want to do this for the "All gear" scenario...
if(unsafe){
  port.foots.stack.noall<- port.foots.stack[[-which(grepl(".All.JGS.PROPORTION", names(port.foots.stack)))]]
} else {
  port.foots.stack.noall<- port.foots.stack[[-which(grepl(".All.JGS.SAFE.PROPORTION", names(port.foots.stack)))]]
}

# Empty stack for results
stack.maxd<- stack()

# Loop
for(i in 1:nlayers(port.foots.stack.noall)){
  lay.use<- port.foots.stack.noall[[i]]
  
  if(unsafe){
    port.name.use<- str_replace_all(names(lay.use), c(".Pot.Trap.JGS.PROPORTION" = "", ".Other.JGS.PROPORTION" = "", ".Gillnet.JGS.PROPORTION" = "", ".Trawl.JGS.PROPORTION" = "", ".Dredge.JGS.PROPORTION" = "", ".Purse.Seine.JGS.PROPORTION" = "", ".Longline.JGS.PROPORTION" = ""))
  } else {
    port.name.use<- str_replace_all(names(lay.use), c(".Pot.Trap.JGS.SAFE.PROPORTION" = "", ".Other.JGS.SAFE.PROPORTION" = "", ".Gillnet.JGS.SAFE.PROPORTION" = "", ".Trawl.JGS.SAFE.PROPORTION" = "", ".Dredge.JGS.SAFE.PROPORTION" = "", ".Purse.Seine.JGS.SAFE.PROPORTION" = "", ".Longline.JGS.SAFE.PROPORTION" = ""))
  }
  pt.use<- unique(ports.geo[str_replace_all(ports.geo$JGSCommunity, "[^[:alnum:]]", "") == str_replace_all(port.name.use, "[^[:alnum:]]", ""),])
  coordinates(pt.use)<- ~Long+Lat
  proj4string(pt.use)<- proj.wgs84
  
  # Get distance from port, then get maximum distance given fishing for max dist and maxdist * 1.5
  dist.rast<- distanceFromPoints(lay.use, pt.use)
  max.d<- dist.rast
  max.d.maxval<- maxValue(raster::mask(dist.rast, lay.use))
  max.d[max.d >= max.d.maxval]<- NA
  max.d<- raster::mask(max.d, lay.use, inverse = T)
  
  max.d.fac<- dist.rast
  max.d.fac[max.d.fac >= max.d.maxval*dist.fac]<- NA
  max.d.fac<- raster::mask(max.d.fac, lay.use, inverse = T)
  
  stack.maxd<- raster::stack(stack.maxd, max.d, max.d.fac)
  print(port.name.use)
}

# Combine
if(unsafe){
  names.stackmaxd<- paste(rep(names(port.foots.stack.noall), each = 2), c("MaxD.JGS.PROPORTION", "1.5xMaxD.JGS.PROPORTION"), sep = "")
} else {
  names.stackmaxd<- paste(rep(names(port.foots.stack.noall), each = 2), c("MaxD.JGS.SAFE.PROPORTION", "1.5xMaxD.JGS.SAFE.PROPORTION"), sep = "")
}
names(stack.maxd)<- names.stackmaxd
names.all<- c(names(port.foots.stack), names.stackmaxd)
check.stack<- raster::stack(port.foots.stack, stack.maxd)
names(check.stack)<- names.all
port.foots.stack<- check.stack

foots.miss

ex.reg<- port.foots.stack[[which(names(port.foots.stack) == paste(foots.miss$Community[[1]], ".", foots.miss$Gear[[1]], ".JGS.PROPORTION", sep =""))]]
ex.max<- port.foots.stack[[which(names(port.foots.stack) == paste(foots.miss$Community[[1]], ".", foots.miss$Gear[[1]], ".JGS.PROPORTIONMaxD.JGS.PROPORTION", sep =""))]]
ex.15max<- port.foots.stack[[which(names(port.foots.stack) == paste(foots.miss$Community[[1]], ".", foots.miss$Gear[[1]], ".JGS.PROPORTION1.5xMaxD.JGS.PROPORTION", sep =""))]]
par(mfrow = c(1, 3))
plot(ex.reg)
plot(ex.max)
plot(ex.15max)

port.foots.stack2<- raster::stack(ex.reg, ex.max, ex.15max)

## Those seem fine???

# Atlantic croaker projections
out.dir<- "~/GitHub/COCA/Results/NormalVoting_BiomassIncPresNoExposure_09222018/"
results<- read_rds(paste(out.dir, "SDMPredictions.rds", sep = "")) %>%
  filter(., COMNAME == "ATLANTIC CROAKER")

# Cold diff and then cold percent difference...
data<- results %>%
  filter(., Proj.Class == "Future_cold.combo.b")
df<- data$Projections[[1]]

dat.sp<- data.frame("x" = df$x, "y" = df$y, "z" = df$Projection)
coordinates(dat.sp)<- ~x+y
proj4string(dat.sp)<- proj.wgs84

pts.rast<- rasterize(dat.sp, port.foots.stack[[1]], field = "z", fun = mean, na.rm = TRUE)
pts.rast<- raster::resample(pts.rast, port.foots.stack[[1]])

res.mean<- vector(length = raster::nlayers(port.foots.stack), mode = "double")
res.mean[]<- NA
res.sd<- vector(length = raster::nlayers(port.foots.stack), mode = "double")
res.sd[]<- NA

i = 2

for(i in 1:raster::nlayers(port.foots.stack2)) {
  lay.use<- port.foots.stack2[[i]]
  
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

# Cold Diff Regular: All zero, Cold Perc Diff (0/0 = NA). Why does this happen? Happens if you have 0 change and a start value of 0 -- yields NaN.
