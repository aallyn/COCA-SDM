library(tidyverse)
library(maptools)
library(raster)
library(rgeos)
library(geosphere)
library(zoo)

suppressWarnings(sapply(list.files(pattern = "[.]R$", path = "~/Dropbox/Andrew/Work/GMRI/AllRFunctions/", full.names = TRUE), source))

climate.dir<- "~/Dropbox/Andrew/Work/GMRI/COCA/Data/ClimateData/climate.sst.proj.grd"
oisst.dir<- "~/Dropbox/Andrew/Work/GMRI/COCA/Data/1981to2015_OISST_Projected.grd"
sp.in<- "~/Dropbox/Andrew/Work/GMRI/AllGIS/"
dates.baseline<- c("2011-10-16", "2012-10-16", "2013-10-16", "2014-10-16", "2015-10-16")
dates.future<- c("2004-10-16", "2012-10-16", "2025-10-16", "2040-10-16", "2055-10-16")
#dates<- c("2015-10-16")
seasonal.mu<- TRUE
season<- "Fall"
model.dat<- "./Data/dat.fitted.SST_01172018.rds"
plot<- TRUE

## Projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

##### Start
## Baseline SSTs
# Empty stack
pred.rast.stack<- stack()

# Add OISST
name.ind<- nlayers(pred.rast.stack)+1
stack0<- stack(oisst.dir)

# Move to monthly?
oisst.min<- gsub("X", "", min(names(stack0)))
oisst.min.date<- as.Date(gsub("[.]", "-", oisst.min))
oisst.max<- gsub("X", "", max(names(stack0)))
oisst.max.date<- as.Date(gsub("[.]", "-", oisst.max))

# Calculate monthly mean temperature -- this would be compared to the sstclim data (monthly climate ensemble)
oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
oisst.dat<- setZ(stack0, oisst.dates)

# Aggregate daily to monthly data
oisst.monthly <- zApply(oisst.dat, by = as.yearmon, mean)

# Climate
name.ind<- nlayers(pred.rast.stack)+1
stack0<- stack(climate.dir)
stack0<- resample(stack0, oisst.monthly[[1]])
clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2099-12-16"), by = "month")
clim.stack<- setZ(stack0, clim.dates)
clim.stack.zind<- getZ(clim.stack)

#Baseline stack, store seasonal means and then average them all
sst.stack<- stack()

years<- c("2011", "2012", "2013", "2014", "2015")

for(i in seq_along(years)){
  
  dates.use<- switch(season,
                     "Fall" = paste(c("Sep", "Oct", "Nov"), rep(years[i]), sep = "."),
                     "Spring" = paste(c("Mar", "Apr", "May"), rep(years[i]), sep = "."))
  
  sst.temp<- calc(oisst.monthly[[which(names(oisst.monthly) %in% dates.use)]], mean)
  
  sst.stack<- stack(sst.stack, sst.temp)
  print(years[i])
}

names(sst.stack)<- paste(season, years, sep = ".")

sst.basemeans<- calc(sst.stack[[c(1,2,3,4,5)]], mean)

pred.rast.stack<- stack(sst.basemeans, sst.stack)
names(pred.rast.stack)[c(1)]<- c("Baseline")

# Climate


# Other predictors
# Add depth
neshelf.bathy<- raster(paste(sp.in, "NEShelf_etopo1_bathy_reclass.tif", sep = ""))
neshelf.bathy[]<- abs(neshelf.bathy[])
proj4string(neshelf.bathy)<- proj.wgs84
depth.temp<- projectRaster(neshelf.bathy, crs = proj.utm)
DEPTH<- resample(depth.temp, pred.rast.stack[[1]])
names(DEPTH)<- "DEPTH"
pred.rast.stack<- stack(pred.rast.stack, DEPTH)

# Shelf_pos
start.pt<- data.frame("x" = -75, "y" = 35)
coordinates(start.pt)<- ~x+y
proj4string(start.pt)<- proj.wgs84
start.pt.utm<- spTransform(start.pt, proj.utm)
shelf_pos<- DEPTH
shelf_pos[]<- NA
names(shelf_pos)<- "SHELF_POS)"

shelf.coords<- data.frame(coordinates(shelf_pos))
coordinates(shelf.coords)<- ~x+y
proj4string(shelf.coords)<- proj.utm
shelf_pos[]<- gDistance(shelf.coords, start.pt.utm, byid = T)/1000
SHELF_POS<- shelf_pos
pred.rast.stack<- stack(pred.rast.stack, SHELF_POS)

# Get these variables
#Mask out points outside of NELME
nelme.rast<- pred.rast.stack[[1]]
nelme.rast[]<- NA
nelme<- readShapePoly(paste("~/Dropbox/Andrew/Work/GMRI/AllGIS/nelme.shp", sep = ""))
proj4string(nelme)<- proj.wgs84
nelme.utm<- spTransform(nelme, proj.utm)
nelme.buff<- gBuffer(nelme.utm, width = 40000)
nelme.rast<- rasterize(nelme.buff, nelme.rast)
pred.rast.stack.m<- mask(pred.rast.stack, mask = nelme.rast, inverse = FALSE)

pred.rast.stack.out<- projectRaster(pred.rast.stack.m, crs = proj.wgs84)


