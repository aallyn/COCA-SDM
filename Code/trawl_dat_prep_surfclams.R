###### Merging SDM and NEVA: Atlantic surfclam
################
######

library(tidyverse)
library(lubridate)
library(sp)
library(sf)
library(rgeos)
library(raster)
library(geosphere)
sst.nc<- "EC_sst_1981_2015_OISST-V2-AVHRR_agg_combined.nc"

## Read in data and join to assessment data
clam.dat<- read.csv("~/GitHub/COCA/Data/clam.model.dat.csv")
names(clam.dat)[23]<- "SPECIES"
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")

clam.dat<- clam.dat %>%
  left_join(., fish.spp, by = "SPECIES") %>%
  drop_na(COMNAME)

# Duplicates?
clam.dat$ID<- paste(clam.dat$cruise, clam.dat$station, clam.dat$tow, sep = ".")
t1<- clam.dat[!duplicated(clam.dat[c("ID", "SVSPP", "catchsex")]),]

# First, we create a null dataset for merging with all species and all trawls
null.dat<- expand.grid("ID" = unique(t1$ID), "SVSPP" = unique(t1$SVSPP))

# Next, create a dataset to match the null but that has species - biomass - abundance for tows where a species was caught.
dat.agg.r<- t1 %>%
  group_by(., ID, SVSPP) %>%
  summarise(., BIOMASS = sum(catch_weight..kg., na.rm = T), ABUNDANCE = sum(catch_number, na.rm = T))

# Merge these together to add in null dataset (absences) for tows where a species was not caught
dat.agg.f<- full_join(dat.agg.r, null.dat)

# Convert NAs to 0, these are where both biomass and abundance are NA
dat.agg.f$bio.abund.na<- ifelse(is.na(dat.agg.f$ABUNDANCE) & is.na(dat.agg.f$BIOMASS), NA, "GOOD")

# Okay, we can make BIOMASS and ABUNDANCE zero where both are NA
dat.agg.f$BIOMASS[is.na(dat.agg.f$bio.abund.na)]<- 0
dat.agg.f$ABUNDANCE[is.na(dat.agg.f$bio.abund.na)]<- 0

# Finally, merge back in the observations with the tow characteristics dataset
tow.dat<- t1[,c(29, 1:21)]
tow.dat<- tow.dat[!duplicated(tow.dat$ID),] # Remove duplicate rows

# Merge tow.dat with dat.agg based on "ID"
dat.full<- full_join(tow.dat, dat.agg.f, by = "ID")
dat.full<- dat.full[order(dat.full$begin_gmt_towdate),]

# For the extractions, going to use projected coordinates UTM zone 19N
# Set projection information
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

# Convert dat.full to spatial points dataframe and we really only need the unique tows
t1<- dat.full[!duplicated(dat.full[c("ID")]),] 

# Remove some of the unneded columns
keep.vec<-  c("ID", "cruise", "vessel", "year", "season", "station", "stratum", "tow", "shg", "begin_gmt_towdate", "end_gmt_towdate", "decdeg_beglat", "decdeg_beglon")
trawl.pred.vals<- dplyr::select(t1, one_of(keep.vec))

# Create spatial points object
trawl.pred.vals.sp<- trawl.pred.vals
coordinates(trawl.pred.vals.sp)<- ~decdeg_beglon+decdeg_beglat
proj4string(trawl.pred.vals.sp)<- proj.wgs84

# Where are these??
nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
st_crs(nelme)<- "+init=epsg:4326"
nelme.sp<- as(nelme, "Spatial")
nelme.sp<- spTransform(nelme.sp, proj.wgs84)

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

plot(nelme.sp)
plot(trawl.pred.vals.sp, col = "red", add = T)

trawl.pred.vals.sp<- spTransform(trawl.pred.vals.sp, proj.utm)

# Getting satellite SST values
dailymu.stack<- raster::stack(paste("./Data/", sst.nc, sep = ""))
dailymu.stack <- raster::rotate(dailymu.stack) # Conversion from 0-360 to -180-180
dailymu.stack<- projectRaster(dailymu.stack, crs = proj.utm)

# Reduce trawl dataset to only those observations that are within the range of SST data
trawl.pred.vals.sp$DATE<- as.Date(format(as.Date(as.POSIXct(trawl.pred.vals.sp$begin_gmt_towdate, format="%m/%d/%y %H:%M", tz="GMT")), "%Y-%m-%d"))
trawl.pred.vals.sp<- trawl.pred.vals.sp[order(trawl.pred.vals.sp$DATE),] 
trawl.pred.vals.sp<- trawl.pred.vals.sp[!is.na(trawl.pred.vals.sp$DATE),]
trawl.pred.vals.sub<- data.frame(trawl.pred.vals.sp[trawl.pred.vals.sp$DATE >= as.Date("1982-01-01") & trawl.pred.vals.sp$DATE <= as.Date("2015-12-31"), ])

# Okay, now with the OISST, can we get the seasonal means for each year...then match years/seasons and extract?
years<- seq(from = 1982, to = 2015, by = 1)
seasonalmu.stack<- stack()

for(i in seq_along(years)){
  summers<- seq(from = as.Date(paste(years[i], "06", "01", sep = "-")), to = as.Date(paste(years[i], "09", "30", sep = "-")), by = 1)
  springs<- seq(from = as.Date(paste(years[i], "05", "01", sep = "-")), to = as.Date(paste(years[i], "06", "30", sep = "-")), by = 1)
  summers.mu<- calc(dailymu.stack[[which(gsub("[.]", "-", gsub("X", "", names(dailymu.stack))) %in% as.character(summers))]], mean)
  spring.mu<- calc(dailymu.stack[[which(gsub("[.]", "-", gsub("X", "", names(dailymu.stack))) %in% as.character(springs))]], mean)
  seasonalmu.stack<- stack(seasonalmu.stack, summers.mu, spring.mu)
  print(years[i])
}

names(seasonalmu.stack)<- paste(c("SUMMER", "SPRING"), rep(years, each = 2), sep = ".")

# Get seasonal temp for each observation
i.temp<- match(paste(toupper(trawl.pred.vals.sub$season), format(as.Date(trawl.pred.vals.sub$DATE), "%Y"), sep = "."), names(seasonalmu.stack)) 
trawl.pred.vals.sub.locations<- trawl.pred.vals.sub[, c("decdeg_beglon", "decdeg_beglat")]
dat.temp<- raster::extract(seasonalmu.stack, trawl.pred.vals.sub.locations) # ~12 minutes

# Bind the i.final vector to our reduced dataset
i.mat <- cbind(1:nrow(trawl.pred.vals.sub), i.temp) 

# Bind result to our sub dataset
trawl.pred.vals.sub$SEASONALMU.OISST<- dat.temp[i.mat]

# Get back to full dataset
temp.pre1982<- data.frame(trawl.pred.vals.sp[trawl.pred.vals.sp$DATE < as.Date("1982-01-01") | trawl.pred.vals.sp$DATE > as.Date("2015-12-31"), ])
temp.pre1982$SEASONALMU.OISST<- NA
trawl.pred.vals<- rbind(temp.pre1982, trawl.pred.vals.sub)

# Getting depth values
neshelf.bathy<- raster("./Data/NEShelf_etopo1_bathy.tiff") 
proj4string(neshelf.bathy)<- proj.wgs84
neshelf.bathy.utm<- projectRaster(neshelf.bathy, crs = proj.utm)

# Use raster to extract raster values at points
trawl.pred.vals$DEPTH<- raster::extract(neshelf.bathy.utm, trawl.pred.vals.sp, method = "bilinear")

# Calculate along shelf position
trawl.pred.vals$SHELF_POS<- distCosine(spTransform(trawl.pred.vals.sp, proj.wgs84), cbind(-75, 35), r=6378137)/1000

# Get ID and added predictor variables
trawl.pred.vals.red<- trawl.pred.vals[,c(1, 12, 16:18)]

# Merge back with trawl catch info (date, time, lat, long, etc)
trawl.dat.full<- full_join(dat.full, trawl.pred.vals.red, by = c("ID"))

# One more modification to get PRESENCE and BIOMASS.MOD data ready
trawl.dat.full$PRESENCE<- ifelse(trawl.dat.full$BIOMASS > 0, 1, trawl.dat.full$BIOMASS) # Create presence/absence vector based on biomass
trawl.dat.full$BIOMASS.MOD<- log(trawl.dat.full$BIOMASS+1)


trawl.dat.full.reduced<- trawl.dat.full[!is.na(trawl.dat.full$SEASONALMU.OISST),] # Keep only data with spring and fall seasonal mean

# Write out model adjusted file -----------------------------------------------------------
saveRDS(trawl.dat.full.reduced, file = "./Data/model.dat.clam.rds")
