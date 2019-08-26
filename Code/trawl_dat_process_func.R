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

trawl_dat_proc<- function(survdat.file = "~/Volumes/Shared/Research/Mills Lab/NMFS Trawl Data/Survdat_Nye_allseason.RData", out.path = "./Data/") {
  ## Details
  # This function reads in an .RData file of NOAA NEFSC bottom trawl survey data and then makes some modifications, including: cutting observations from southeast and Canada strata, adding in "absence" observations, summarizing biomass/abundance across to species at each tow (aggregating across sex and weight and length) and adds a few unique ID comumns. The processed datafile is then saved in out.path as "survdat_processed.dat" with the year and month and day date when the function was executed. 
  
  # Args:
  # survdat.file = Path to survey data file 
  # out.path = Path to save processsed rds data file
  
  # Returns: RDS Data file, which is also saved in folder specified by out.path
  
  ## Start function
  # Preliminaries -----------------------------------------------------------
  # Load libraries, using package_check to download if it is missing
  libraries.needed<- c("tidyverse", "sp", "raster", "geosphere")
  library_check(libraries.needed)
  
  # For debugging
  if(FALSE){
    survdat.file = "./Data/Survdat_Nye2016.RData"
    survdat.file = "/Volumes/Shared/Research/Mills Lab/NMFS Trawl Data/Survdat_Nye_allseason.RData"
    sst = "/Volumes/Shared/Research/Mills Lab/SST/RawData/OISST.grd"
    out.path = "~/GitHub/COCA/Data/"
  }
  
  # Load in the data
  load(survdat.file)
  
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
  
  # Create a new TRAWL.ID column for unique ID
  dat.full$TRAWL.ID<- paste(dat.full$EST_YEAR, dat.full$SEASON, dat.full$STRATUM, dat.full$STATION, sep = ".")
  dat.full<- dat.full[order(dat.full$EST_YEAR, dat.full$EST_MONTH, dat.full$EST_DAY),]

  # One more modification to get PRESENCE and BIOMASS.MOD data ready
  dat.full$PRESENCE<- ifelse(dat.full$BIOMASS > 0, 1, dat.full$BIOMASS) # Create presence/absence vector based on biomass
  dat.full$BIOMASS.MOD<- log(dat.full$BIOMASS+1)
  
  dat.full$DATE<- as.Date(paste(dat.full$EST_YEAR, dat.full$EST_MONTH, dat.full$EST_DAY, sep = "-"))
  
  # Write out model adjusted file -----------------------------------------------------------
  saveRDS(dat.full, file = paste(out.path, "survdat_processed.dat", gsub("-", "", format(Sys.time(), "%Y-%m-%d")), ".rds", sep = ""))
  
  ## End function
}

