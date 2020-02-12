#####
## Climate projected temperatures 
#####


# PRELIMINARIES -----------------------------------------------------------
# Source "library_check" helper function to install/load required libraries
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")
source(paste(lab.func.path, "GeneralHelpers.R", sep = ""))
library_check(c("tidyverse", "here", "raster"))


# PRE-PROCESSING: MEAN, 5TH AND 95H PERCENTILES ---------------------------
## Climate projections file -- has each ensemble member
clim.proj<- readRDS(paste(res.data.path, "CMIP5_SST/ProcessedSSTProjections.rds", sep = ""))

## NELME shapefile
nelme<- st_read(paste(res.data.path, "Shapefiles/NELME_regions/NELME_sf.shp", sep = ""))
proj.wgs84<- "+init=epsg:4326"
st_crs(nelme)<- proj.wgs84

## Calculate the mean, 5th and 95th across the different ensemble members by scenario/year/month/x/y
clim.summs<- clim.proj %>%
  separate(Year, into = c("Year", "Month", "Day")) %>%
  group_by(Scenario, Year, Month, x, y) %>%
  summarize("Mean" = mean(SST, na.rm = TRUE),
            "Pct5th" = quantile(SST, probs = c(0.05), na.rm = TRUE, names = FALSE),
            "Pct95th" = quantile(SST, probs = c(0.95), na.rm = TRUE, names = FALSE))

## Nest the data to make things a bit easier to work with
clim.summs<- clim.summs %>%
  group_by(Scenario, Year, Month) %>%
  nest()

## Function to convert dataframes to raster layers and use nelme to mask the data
df_to_rast<- function(df, stat, mask = nelme) {
  if(FALSE){
    df<- clim.summs$data[[1]]
    stat<- "Mean.SST"
  }
  
  df.temp<- df %>%
    dplyr::select(x, y, stat)
  
  rast.temp<- rasterFromXYZ(df.temp)
  rast.out<- mask(rast.temp, nelme)
  return(rast.out)
}

## Apply the function to each of the nested dataframes
clim.summs<- clim.summs %>%
  mutate(., "RasterStack.Mean" = map2(data, "Mean", df_to_rast),
         "RasterStack.Pct05" = map2(data, "Pct5th", df_to_rast),
         "RasterStack.Pct95" = map2(data, "Pct95th", df_to_rast))

## Now, need to create a stack for the means, pct05, pct95
# 8.5 First
rcp85<- clim.summs[clim.summs$Scenario == "RCP85", ]
rcp85.mu<- raster::stack(rcp85$RasterStack.Mean)
names(rcp85.mu)<- paste(rcp85$Month, rcp85$Year, sep = "-")
proj4string(rcp85.mu)<- proj.wgs84
writeRaster(rcp85.mu, paste(res.data.path, "CMIP5_SST/RCP85_mu.grd", sep = ""))

rcp85.pct05<- raster::stack(rcp85$RasterStack.Pct05)
names(rcp85.pct05)<- paste(rcp85$Month, rcp85$Year, sep = "-")
proj4string(rcp85.pct05)<- proj.wgs84
writeRaster(rcp85.pct05, paste(res.data.path, "CMIP5_SST/RCP85_5th.grd", sep = ""))

rcp85.pct95<- raster::stack(rcp85$RasterStack.Pct95)
names(rcp85.pct95)<- paste(rcp85$Month, rcp85$Year, sep = "-")
proj4string(rcp85.pct95)<- proj.wgs84
writeRaster(rcp85.pct95, paste(res.data.path, "CMIP5_SST/RCP85_95th.grd", sep = ""))

# Now, 4.5
rcp45<- clim.summs[clim.summs$Scenario == "RCP45", ]
rcp45.mu<- raster::stack(rcp45$RasterStack.Mean)
names(rcp45.mu)<- paste(rcp45$Month, rcp45$Year, sep = "-")
proj4string(rcp45.mu)<- proj.wgs84
writeRaster(rcp45.mu, paste(res.data.path, "CMIP5_SST/RCP45_mu.grd", sep = ""))

rcp45.pct05<- raster::stack(rcp45$RasterStack.Pct05)
names(rcp45.pct05)<- paste(rcp45$Month, rcp45$Year, sep = "-")
proj4string(rcp45.pct05)<- proj.wgs84
writeRaster(rcp45.pct05, paste(res.data.path, "CMIP5_SST/RCP45_5th.grd", sep = ""))

rcp45.pct95<- raster::stack(rcp45$RasterStack.Pct95)
names(rcp45.pct95)<- paste(rcp45$Month, rcp45$Year, sep = "-")
proj4string(rcp45.pct95)<- proj.wgs84
writeRaster(rcp45.pct95, paste(res.data.path, "CMIP5_SST/RCP45_95th.grd", sep = ""))