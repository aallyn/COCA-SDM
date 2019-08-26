#####
## FVCOM Bottom Temperatures Code
#####

# Set path to "fvcom_extract_func.R" and source the functions included there (libraries_check and fvcom_extract)
func.path<- "~/GitHub/COCA/"
source(paste(func.path, "fvcom_extract_func.R", sep = ""))

# Run the main function -- set up right now to get monthly mean hindcasted bottom temperatures for Gulf of Maine
temps.out<- fvcom_extract(threddsURL = "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean", lonW = -82, lonE = -61.75, latS = 24, latN = 46.25, spacing = 0.2, date1 = NULL, date2 = NULL, variable = "temp", out.path = func.path)

# Plot one
plot(temps.out[[1]])

# Read back in the raster stack....
stack.path<- paste(func.path, "FVCOMtemp.grd", sep = "")
temps.stack<- raster::stack(stack.path)

dev.off()
plot(temps.stack[[1]])
