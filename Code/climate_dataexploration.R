###############
### COCA Climate Data Exploration ###
###############

## Libraries
library(ncdf4)
library(raster)
library(maptools)
library(ggplot2)
library(tidyverse)
library(rgeos)
library(zoo)
library(viridis)

## Some spatial stuff for data visualiztion
# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

#Bounds
xlim.use<- c(-77, -65)
ylim.use<- c(35.05, 45.2)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us <- raster::getData("GADM",country="USA",level=1)
us.states <- us[us$NAME_1 %in% states,]
us.states <- gSimplify(us.states, tol = 0.075, topologyPreserve = TRUE)
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]
ca.provinces <- gSimplify(ca.provinces, tol = 0.075, topologyPreserve = TRUE)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

#########
#### Applying anomalies -- mean
#########
## Inspect the file using the ncdf4 libaray
## Get sst anomaly 
sst.anom.temp<- raster::stack("~/GitHub/COCA/Data/SST.CMIP5.1982-2099.anom.nc", varname = "sstanom")
sst.anom<- raster::rotate(sst.anom.temp)

## Get oisst data
oisst.dat.temp<- raster::stack("~/Dropbox/Andrew/Work/GMRI/Projects/AllData/EC_sst_1981_2015_OISST-V2-AVHRR_agg_combined.nc")
oisst.dat<- raster::rotate(oisst.dat.temp)

# Need to get climatology from the OISST data -- set up OISST stack as time series
oisst.min<- gsub("X", "", min(names(oisst.dat)))
oisst.min.date<- as.Date(gsub("[.]", "-", oisst.min))
oisst.max<- gsub("X", "", max(names(oisst.dat)))
oisst.max.date<- as.Date(gsub("[.]", "-", oisst.max))

# Calculate monthly mean temperature -- this would be compared to the sstclim data (monthly climate ensemble)
oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
oisst.dat<- setZ(oisst.dat, oisst.dates)

# Aggregate daily to monthly data
oisst.monthly <- zApply(oisst.dat, by = as.yearmon, mean)

# Lets check that
test.dat<- subset(oisst.dat, which(getZ(oisst.dat) >="1981-09-01" & getZ(oisst.dat) <= "1981-09-30"))
sept1981.mu<- calc(test.dat, mean)
plot(oisst.monthly[[1]]-sept1981.mu)

# Everything seems fine there, now need the monthly average for each month across baseline years (1982-2011)
dates<- getZ(oisst.monthly)
subset.vec<- which(dates > "Dec 1981" & dates < "Jan 2012", arr.ind = TRUE)
oisst.monthly.sub<- oisst.monthly[[subset.vec]]
oisst.monthly.sub<- setZ(oisst.monthly.sub, dates[subset.vec])

oisst.clim<- stack()
months<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

for(i in seq_along(months)) {
  # Get all the monthly raster
  stack.temp<- subset(oisst.monthly.sub, which(grepl(months[i], names(oisst.monthly.sub))))
  month.clim<- calc(stack.temp, fun = mean)
  oisst.clim<- stack(oisst.clim, month.clim)
  names(oisst.clim)[i]<- months[i]
}

# Check that
test.dat<- subset(oisst.monthly.sub, which(grepl("Jan", names(oisst.monthly.sub))))
jan.mu<- calc(test.dat, mean)
summary(oisst.clim[[1]] - jan.mu)

# Looks good -- time to apply the anomalies to the oisst.clim
oisst.clim.coarse<- raster::resample(oisst.clim, sst.anom)
names(oisst.clim.coarse)<- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Okay, now good to apply the anomalies from the climate models to climatology and get raw values
sst.model<- stack()

for(i in 1:nlayers(sst.anom)) {
  index.match<- which(gsub("X", "", names(oisst.clim.coarse)) == unlist(strsplit(names(sst.anom)[i], "[.]"))[2])
  rast.temp<- oisst.clim.coarse[[index.match]] + sst.anom[[i]]
  sst.model<- stack(sst.model, rast.temp)
  names(sst.model)[i]<- names(sst.anom)[i]
}

# One more step, need to fill in the missing coastline raster cell values.... Function below, i is corresponding to a three by three window
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

# Now apply that function to each raster stack
for(i in 1:nlayers(sst.model)) {
  new.rast<- focal(sst.model[[i]], w = matrix(1, 3, 3), fun = fill.na, pad = TRUE, na.rm = FALSE)
  sst.model[[i]]<- new.rast
}

sst.model.proj<- projectRaster(sst.model, crs = proj.utm)
names(sst.model.proj)<- names(sst.anom)
sst.model<- projectRaster(sst.model.proj, crs = proj.wgs84)
writeRaster(sst.model.proj, filename = paste(data.dir, "climate.sst.proj.grd", sep = ""), format = "raster", overwrite = TRUE)

# Looking at seasonal changes
fall<- c("X2011.09.16", "X2011.10.16", "X2011.11.16", "X2012.09.16", "X2012.10.16", "X2012.11.16", "X2013.09.16", "X2013.10.16", "X2013.11.16", "X2014.09.16", "X2014.10.16", "X2014.11.16", "X2015.09.16", "X2015.10.16", "X2015.11.16")
spring<- c("X2011.03.16", "X2011.04.16", "X2011.05.16", "X2012.03.16", "X2012.04.16", "X2012.05.16", "X2013.03.16", "X2013.04.16", "X2013.05.16", "X2014.03.16", "X2014.04.16", "X2014.05.16", "X2015.03.16", "X2015.04.16", "X2015.05.16")

base.fall<- mean(sst.model[[which(names(sst.model) %in% fall)]])
base.spring<- mean(sst.model[[which(names(sst.model) %in% spring)]])

fall.f<- c("X2055.09.16", "X2055.10.16", "X2055.11.16")
fall.diffs<- mean(sst.model[[which(names(sst.model) %in% fall.f)]]) - base.fall
fall.diffs.df<- as.data.frame(fall.diffs, xy = T)
fall.diffs.df$Season<- rep("Fall", nrow(fall.diffs.df))

spring.f<- c("X2055.03.16", "X2055.04.16", "X2055.05.16")
spring.diffs<- mean(sst.model[[which(names(sst.model) %in% spring.f)]]) - base.spring
spring.diffs.df<- as.data.frame(spring.diffs, xy = T)
spring.diffs.df$Season<- rep("Spring", nrow(spring.diffs.df))

diffs.df<- bind_rows(fall.diffs.df, spring.diffs.df)

ggplot() + 
  geom_tile(data = diffs.df, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(name = "SST", option = "viridis", na.value = "white") +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black")) +
  facet_wrap(~Season)

# Line plot
diffs.neslme<- diffs.df


#####
## Min and Max
#####
## Get sst anomaly 
sst.anom.temp<- raster::stack("~/GitHub/COCA/Data/SST.CMIP5.1982-2099.anom.nc", varname = "sstpct05")
sst.anom<- raster::rotate(sst.anom.temp)

# Okay, now good to apply the anomalies from the climate models to climatology and get raw values
sst.model<- stack()

for(i in 1:nlayers(sst.anom)) {
  index.match<- which(gsub("X", "", names(oisst.clim.coarse)) == unlist(strsplit(names(sst.anom)[i], "[.]"))[2])
  rast.temp<- oisst.clim.coarse[[index.match]] + sst.anom[[i]]
  sst.model<- stack(sst.model, rast.temp)
  names(sst.model)[i]<- names(sst.anom)[i]
}

# One more step, need to fill in the missing coastline raster cell values.... Function below, i is corresponding to a three by three window
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

# Now apply that function to each raster stack
for(i in 1:nlayers(sst.model)) {
  new.rast<- focal(sst.model[[i]], w = matrix(1, 3, 3), fun = fill.na, pad = TRUE, na.rm = FALSE)
  sst.model[[i]]<- new.rast
}

sst.model.proj<- projectRaster(sst.model, crs = proj.utm)
names(sst.model.proj)<- names(sst.anom)
sst.model<- projectRaster(sst.model.proj, crs = proj.wgs84)
writeRaster(sst.model.proj, filename = paste("~/GitHub/COCA/Data/", "climate.sst.proj.pct05.grd", sep = ""), format = "raster", overwrite = TRUE)

## Get sst anomaly 
sst.anom.temp<- raster::stack("~/GitHub/COCA/Data/SST.CMIP5.1982-2099.anom.nc", varname = "sstpct95")
sst.anom<- raster::rotate(sst.anom.temp)

# Okay, now good to apply the anomalies from the climate models to climatology and get raw values
sst.model<- stack()

for(i in 1:nlayers(sst.anom)) {
  index.match<- which(gsub("X", "", names(oisst.clim.coarse)) == unlist(strsplit(names(sst.anom)[i], "[.]"))[2])
  rast.temp<- oisst.clim.coarse[[index.match]] + sst.anom[[i]]
  sst.model<- stack(sst.model, rast.temp)
  names(sst.model)[i]<- names(sst.anom)[i]
}

# One more step, need to fill in the missing coastline raster cell values.... Function below, i is corresponding to a three by three window
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

# Now apply that function to each raster stack
for(i in 1:nlayers(sst.model)) {
  new.rast<- focal(sst.model[[i]], w = matrix(1, 3, 3), fun = fill.na, pad = TRUE, na.rm = FALSE)
  sst.model[[i]]<- new.rast
}

sst.model.proj<- projectRaster(sst.model, crs = proj.utm)
names(sst.model.proj)<- names(sst.anom)
sst.model<- projectRaster(sst.model.proj, crs = proj.wgs84)
writeRaster(sst.model.proj, filename = paste("~/GitHub/COCA/Data/", "climate.sst.proj.pct95.grd", sep = ""), format = "raster", overwrite = TRUE)



