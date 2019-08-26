#####
## Functions for processing NOAA NEFSC bottom trawl survey data 
#####

# Helper Function ---------------------------------------------------------
library_check<- function(libraries) {
  ## Details
  # This function will check for and then either load or install any libraries needed to run subsequent functions
  
  # Args:
  # libraries = Vector of required library names
  
  # Returns: NA, just downloads/loads libraries into current work space.
  
  ## Start function
  lapply(libraries, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
  ## End function
}

trawl_dat_prep<- function(survdat.file = "~/Volumes/Shared/Research/MillsLab/NMFS Trawl Data/Survdat_Nye_allseason.RData", sst = "~/Volumes/Shared/Research/MillsLab/SST/RawData/OISST.grd", out.path = "./Data/") {
  ## Details
  # This function reads in an .RData file of NOAA NEFSC bottom trawl survey data and then makes some modifications to produce a datafile ready for fitting with SDM, such that for every tow there is a record of every species (presence, log(biomass+1) given presence) and seasonal SST average and bottom depth. Seasonal SST averages are calculated from the supplied .nc OISST time series file. After cleaning the raw data and collecting the predictor variables, the function saves the model ready file as "model.dat.rds."
  
  # Args:
  # survdat.file = Path to survey data file 
  # sst = Path to sea surface temperature data, hosted on the shared drive
  # out.path = Path to save processsed rds data file
  
  # Returns: RDS Data file, which is also saved in folder specified by out.path
  
  ## Start function
  # Preliminaries -----------------------------------------------------------
  # Load libraries, using package_check to download if it is missing
  libraries.needed<- c("tidyverse", "sp", "raster", "geosphere")
  library_check(libraries.needed)
  
  # Load in the data
  load(paste("./Data/", survdat.file, sep = ""))
  
  # Data cleaning -----------------------------------------------------------
  # Remove tows from Canadian waters and SEFSC 
  strata.ca<- c(01350, 1351, 1310, 1320, 1410, 1420, 1490, 1990, 1410, 1420, 1490, 5440, 5480, 5430) # Canada
  strata.se<- unique(survdat$STRATUM[survdat$STRATUM > 3990]) #South east strata and other survey strata
  strata.delete<- c(strata.ca, strata.se)
  t1<- survdat[!survdat$STRATUM %in% strata.delete,]
  
  # There are duplicated rows for tow, species, biomass, abundance, catchsex for different lengths of individuals. We are only interested in variability in biomass for this work. So, we can reduce the dataset by ignoring these repeating rows.
  t2<- t1[!duplicated(t1[c("ID", "SVSPP", "CATCHSEX")]),]
  
  # Now, we have unique tow-species-biomass-abundance records for tows were a species was caught. However, we also need to abesences, so we need biomass and abundance of zero for each species when they were not caught at a specific tow.
  # First, we create a null dataset for merging with all species and all trawls
  null.dat<- expand.grid("ID" = unique(t2$ID), "SVSPP" = unique(t2$SVSPP))
  
  # Next, create a dataset to match the null but that has species - biomass - abundance for tows where a species was caught.
  dat.agg.r<- t2 %>%
    group_by(., ID, SVSPP) %>%
    summarise(., BIOMASS = sum(BIOMASS, na.rm = T), ABUNDANCE = sum(ABUNDANCE, na.rm = T))
  
  # Merge these together to add in null dataset (absences) for tows where a species was not caught
  dat.agg.f<- full_join(dat.agg.r, null.dat)
  
  # Convert NAs to 0, these are where both biomass and abundance are NA
  dat.agg.f$bio.abund.na<- ifelse(is.na(dat.agg.f$ABUNDANCE) & is.na(dat.agg.f$BIOMASS), NA, "GOOD")
  
  # Okay, we can make BIOMASS and ABUNDANCE zero where both are NA
  dat.agg.f$BIOMASS[is.na(dat.agg.f$bio.abund.na)]<- 0
  dat.agg.f$ABUNDANCE[is.na(dat.agg.f$bio.abund.na)]<- 0
  
  # Finally, merge back in the observations with the tow characteristics dataset
  tow.dat<- t2[,c(1:29, 35:43)]
  tow.dat<- tow.dat[!duplicated(tow.dat$ID),] # Remove duplicate rows
  
  # Merge tow.dat with dat.agg based on "ID"
  dat.full<- full_join(tow.dat, dat.agg.f, by = c("ID"))
  
  # Extracting environmental variables at observation locations -----------------------------------------------------------
  # Create a new TRAWL.ID column for unique ID
  dat.full$TRAWL.ID<- paste(dat.full$EST_YEAR, dat.full$SEASON, dat.full$STRATUM, dat.full$STATION, sep = ".")
  dat.full<- dat.full[order(dat.full$EST_YEAR, dat.full$EST_MONTH, dat.full$EST_DAY),]
  
  # For the extractions, going to use projected coordinates UTM zone 19N
  # Set projection information
  proj.wgs84<- CRS("+init=epsg:4326") #WGS84
  proj.utm<- CRS("+init=epsg:2960") #UTM 19
  
  # Convert dat.full to spatial points dataframe and we really only need the unique tows
  t1<- dat.full[!duplicated(dat.full[c("ID")]),] 
  
  # Remove some of the unneded columns
  keep.vec<-  c("ID", "STATUS_CODE", "EST_YEAR", "EST_MONTH", "EST_DAY", "EST_JULIAN_DAY", "EST_TIME", "SEASON", "SVVESSEL", "TOWDUR", "RPM", "BOTSPEED", "BOTTEMP", "SURFTEMP", "BOTSALIN", "SURFSALIN", "STRATUM", "DECDEG_BEGLAT", "DECDEG_BEGLON", "TRAWL.ID")
  trawl.pred.vals<- dplyr::select(t1, one_of(keep.vec))
  
  # Create spatial points object
  trawl.pred.vals.sp<- trawl.pred.vals
  coordinates(trawl.pred.vals.sp)<- ~DECDEG_BEGLON+DECDEG_BEGLAT
  proj4string(trawl.pred.vals.sp)<- proj.wgs84
  trawl.pred.vals.sp<- spTransform(trawl.pred.vals.sp, proj.utm)
  
  # Getting satellite SST values
  dailymu.stack<- raster::stack(paste("./Data/", sst.nc, sep = ""))
  dailymu.stack <- raster::rotate(dailymu.stack) # Conversion from 0-360 to -180-180
  dailymu.stack<- projectRaster(dailymu.stack, crs = proj.utm)
  
  # Reduce trawl dataset to only those observations that are within the range of SST data
  trawl.pred.vals.sp$DATE<- as.Date(paste(trawl.pred.vals.sp$EST_YEAR, trawl.pred.vals.sp$EST_MONTH, trawl.pred.vals.sp$EST_DAY, sep = "-"))
  trawl.pred.vals.sub<- data.frame(trawl.pred.vals.sp[trawl.pred.vals.sp$DATE >= as.Date("1982-01-01"), ])
  
  # Okay, now with the OISST, can we get the seasonal means for each year...then match years/seasons and extract?
  years<- seq(from = 1982, to = 2015, by = 1)
  seasonalmu.stack<- stack()
  
  for(i in seq_along(years)){
    falls<- seq(from = as.Date(paste(years[i], "09", "01", sep = "-")), to = as.Date(paste(years[i], "11", "30", sep = "-")), by = 1)
    springs<- seq(from = as.Date(paste(years[i], "03", "01", sep = "-")), to = as.Date(paste(years[i], "05", "31", sep = "-")), by = 1)
    fall.mu<- calc(dailymu.stack[[which(gsub("[.]", "-", gsub("X", "", names(dailymu.stack))) %in% as.character(falls))]], mean)
    spring.mu<- calc(dailymu.stack[[which(gsub("[.]", "-", gsub("X", "", names(dailymu.stack))) %in% as.character(springs))]], mean)
    seasonalmu.stack<- stack(seasonalmu.stack, fall.mu, spring.mu)
    print(years[i])
  }
  
  names(seasonalmu.stack)<- paste(c("FALL", "SPRING"), rep(years, each = 2), sep = ".")
  
  # Get seasonal temp for each observation
  i.temp <- match(paste(trawl.pred.vals.sub$SEASON, trawl.pred.vals.sub$EST_YEAR, sep = "."), names(seasonalmu.stack)) 
  trawl.pred.vals.sub.locations<- trawl.pred.vals.sub[, c("DECDEG_BEGLON", "DECDEG_BEGLAT")]
  dat.temp <- raster::extract(seasonalmu.stack, trawl.pred.vals.sub.locations) # ~12 minutes
  
  # Bind the i.final vector to our reduced dataset
  i.mat <- cbind(1:nrow(trawl.pred.vals.sub), i.temp) 
  
  # Bind result to our sub dataset
  trawl.pred.vals.sub$SEASONALMU.OISST<- dat.temp[i.mat]
  
  # Get back to full dataset
  temp.pre1982<- data.frame(trawl.pred.vals.sp[trawl.pred.vals.sp$DATE < as.Date("1982-01-01"), ])
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
  trawl.pred.vals.red<- trawl.pred.vals[,c(1, 23:25)]
  
  # Merge back with trawl catch info (date, time, lat, long, etc)
  trawl.dat.full<- full_join(dat.full, trawl.pred.vals.red, by = c("ID"))
  
  # One more modification to get PRESENCE and BIOMASS.MOD data ready
  trawl.dat.full$PRESENCE<- ifelse(trawl.dat.full$BIOMASS > 0, 1, trawl.dat.full$BIOMASS) # Create presence/absence vector based on biomass
  trawl.dat.full$BIOMASS.MOD<- log(trawl.dat.full$BIOMASS+1)
  
  trawl.dat.full$DATE<- as.Date(paste(trawl.dat.full$EST_YEAR, trawl.dat.full$EST_MONTH, trawl.dat.full$EST_DAY, sep = "-"))
  trawl.dat.full.reduced<- trawl.dat.full[trawl.dat.full$DATE >= "1982-01-01",] # Keep only data with spring and fall seasonal mean
  
  # Write out model adjusted file -----------------------------------------------------------
  saveRDS(trawl.dat.full.reduced, file = "./Data/model.dat.rds")
  
  #########
  ## End
  #########
}

