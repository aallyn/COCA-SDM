#####
## FVCOM Data Extraction Functions
#####

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

fvcom_extract<- function(threddsURL = "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean", lonW = -82, lonE = -61.75, latS = 24, latN = 46.25, spacing, date1 = NULL, date2 = NULL, variable, out.path) {
  # This function extracts FVCOM data from the threddsURL. After extracting the data (and subsetting if dates are specified), it saves a the raw data as a large matrix and also as a regular gridded raster stack. To retrieve the full time series takes about ~5 minutes.
  
  # Args:
  # threddsURL = Link to FVCOM THREDDS URL. Default is "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean" A quick note here, this path is actually found at the page "http://www.smast.umassd.edu:8080/thredds/hindcasts.html?dataset=fvcom/hindcasts/30yr_gom3/mean" and under the "Access" options. I think if you wanted to get something else, say the FVCOM forecasts, you would pull the THREDDS URL from this page: http://www.smast.umassd.edu:8080/thredds/forecasts.html?dataset=gom3_nocache. You would then probably need to do some function modifications/checking to make sure everything still works as you expect. 
  # lonW = Left bounding box decimal degree coordinate for raster stack
  # lonE = Right bounding boc decimal degree coordinate for raster stack
  # latS = Bottom bounding box decimal degree coordinate for raster stack
  # latN = Upper bounding box decimal degree coordinate for raster stack
  # spacing = Grid spacing for output raster stack. The FVCOM data is not on a regular grid, which is a requirement of raster stacks. To get around this, the function interpolates the FVCOM data onto a raster layer, defined by lon, lat and spacing. 
  # date1 = Either a date string ("YYYY-MM-DD") to set beginning of time series subset, or NULL to retrieve full time series. Default is NULL.
  # date2 = Either a date string ("YYYY-MM-DD") to set end of time series subset, or NULL to retrieve full time series. Default is NULL.
  # variable = FVCOM variable data to extract. Default is temp, could also be "salinity" -- not sure on others.
  # out.path = Path to save large matrix and raster stack FVCOM data
  
  # Returns: Raster stack of FVCOM variable data. 
  
  # Preliminaries
  library_check(c("ncdf4", "raster", "fields"))
  
  # Debugging
  if(FALSE) {
    threddsURL<- "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean"
    lonW<- -82
    lonE<- -61.75
    latS<- 24
    latN<- 46.25
    spacing<- 0.15
    date1<- NULL
    date2<- NULL
    variable<- "temp"
    out.path<- "~/Desktop/"
  }
  
  ## Start
  # Open file connection, extract lon/lat info
  nc<- nc_open(threddsURL)
  lon.nc<- ncvar_get(nc, varid = "lon")
  nlons<- dim(lon.nc)
  lat.nc<- ncvar_get(nc, varid = "lat")
  nlats<- dim(lat.nc)
  
  # Extract available dates from netCDF file
  ncdates<- nc$dim$time$vals
  ncdates<- as.Date(ncdates,origin = '1858-11-17') #available time points in nc
  
  # Subset full time series?
  if(is.null(date1) && is.null(date2)) {
    date1<- min(ncdates)
    date2<- max(ncdates)
  }
  
  if (class(date1) == 'Date'){
    # Get index of nearest time point
    date1indx = which.min(abs(date1 - ncdates)) 	
  } else if (class(date1) == 'character'){
    # Convert to a Date object first
    date1 = as.Date(date1)
    date1indx = which.min(abs(date1 - ncdates)) 
  }
  if (missing(date2)) {
    # If date2 isn't specified, reuse date1
    date2indx = which.min(abs(date1 - ncdates)) 
    cat('Only 1 date specified\n')
  } else {
    if (class(date2) == 'Date'){
      # If date2 exists, get index of nearest time point to date2
      date2indx = which.min(abs(date2 - ncdates)) 		
    } else if (class(date2) == 'character'){
      date2 = as.Date(date2)
      date2indx = which.min(abs(date2 - ncdates))
    }
  }
  
  # Get number of time steps to extract
  ndates <- (date2indx - date1indx) + 1 #get number of time steps to extract
  
  # Get variable data
  var.fvcom<- ncvar_get(nc, varid = variable, start =  c(1, 45, 1))
  
  # Little clean up to save all data as a large 
  fvcom.df.out<- data.frame("lon" = lon.nc, "lat" = lat.nc, var.fvcom)
  colnames(fvcom.df.out)[3:ncol(fvcom.df.out)]<- paste(variable, ncdates, sep = ".")
  write_csv(fvcom.df.out, paste(out.path, "FVCOM", variable, ".csv", sep = ""))
  
  # Now need to convert to a raster. In the var.fvcom matrix, each row is a point and the columns are the different dates. Need to extract and put into raster stack. One issue with this, though, is that the FVCOM data are not on a perfect grid -- so some interpolation is going to be needed. To do this, set up an empty raster layer
  rast.lon<- seq(from = lonE, to = lonW, by = -spacing)
  rast.lat<- seq(from = latS, to = latN, by = spacing)
  e<- extent(c(lonW, lonE, latS, latN))
  rast.temp<- raster(e, nrows = length(rast.lat), ncols = length(rast.lon))
  
  stack0<- stack()
  
  for (i in 1:ncol(var.fvcom)) {
    t1<- data.frame("long" = lon.nc, "lat" = lat.nc, "z" = var.fvcom[,i])
    rast<- rasterize(t1[,1:2], rast.temp, t1[,3], fun = mean)
    stack0<- stack(stack0, rast)
    print(paste(ncdates[i], " is done", sep = ""))
  }
  
  # Raster stack spatial projection and names
  proj4string(stack0)<- CRS("+init=epsg:4326") #WGS84
  names(stack0)<- as.character(ncdates)
  
  # Save it and return it
  writeRaster(stack0, paste(out.path, "FVCOM", variable, ".grd", sep = ""))
  return(stack0)
}