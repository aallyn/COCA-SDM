####### 
### For Kathy: Filtering by Jon's species
#######

# Helper function to install any needed missing packages
package_check<- function(packages) {
  lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}

# Run it to install "tidyverse" -- may take a bit if it isn't installed already
package_check("tidyverse")

# Set paths to trawl data and assessment fish species -- you will need to update these!!
trawl.dat.path<- "~/GitHub/COCA/Data/Survdat_Nye2016.RData"
hare.spp.path<- "~/GitHub/COCA/Data/Assesmentfishspecies.csv"

# Read in the trawl data and assessment data 
trawl.dat.full<- load(trawl.dat.path)
hare.spp<- read.csv(hare.spp.path)

# Filter trawl data and only include species that were in Jon's assessment (based on SVSPP species code)
trawl.dat.sub<- trawl.dat.full[trawl.dat.full$SVSPP %in% hare.spp$SVSPP,]

# Add Hare columns, which include species common name, to the trawl data
trawl.dat.sub<- left_join(trawl.dat.sub, hare.spp, by = "SVSPP")
