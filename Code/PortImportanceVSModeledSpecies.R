#####
## COCA weighted landings
#####
library(tidyverse)

# Load in model projections and fit
out.dir<- "~/GitHub/COCA/Results/NormalVoting_BiomassIncPresNoExposure_09222018/"
dat.full<- read_rds(paste(out.dir, "SDMPredictions.rds", sep = "")) 

# Model fit stats and filter out models that had less than 0.7 or have low certainty according to expert assessment
mod.results<- read.csv(paste(out.dir, "mod.results.csv", sep = "")) %>%
  dplyr::select(., -X)
dat.full<- dat.full %>%
  left_join(., mod.results, by = c("COMNAME", "SEASON"))
dat.full$AUC.SDM[is.na(dat.full$AUC.SDM)]<- 0

dat.sub<- dat.full %>%
  filter(., AUC.SDM >= 0.70)

dat.mod.results<- dat.sub %>%
  dplyr::select(., COMNAME, SEASON, DevExp.P, DevExp.B, AUC.SDM, Calib.SDM, RMSE.SDM.P, RMSE.NEVA.B) %>%
  group_by(., COMNAME) %>%
  summarize_if(., is.numeric, mean, na.rm = T) %>%
  ungroup() %>%
  mutate(., "Match.Col" = gsub(" ", ".", COMNAME))

# Add model results to landings dataset
landings.dat<- read.csv("~/GitHub/COCA/Data/port_dealer_spp.csv")



%>%
  gather(., "Match.Col", "Landings", -Port.short, -Port.state, -Port.long)
landings.dat$Port.short<- gsub("[.]", "_", landings.dat$Port.short)
landings.dat$Port.long<- paste(landings.dat$Port.short, landings.dat$Port.state, sep = ".")

landings.total<- landings.dat %>%
  group_by(., Port.long) %>%
  summarize_at(., "Landings", sum, na.rm = T)
colnames(landings.total)[2]<- "Total.Landings"

landings.prop<- landings.dat %>%
  left_join(., landings.total, by = "Port.long") %>%
  mutate(., "Prop.Landings" = round(Landings/Total.Landings, 4))

## Species 
landings.prop.modeled.spp<- landings.prop %>%
  left_join(., dat.mod.results, by = "Match.Col") %>%
  drop_na(., AUC.SDM) %>%
  group_by(., Port.long)

landings.prop.modeled.spp$Port.Match<- as.character(landings.prop.modeled.spp$Port.long)
landings.prop.modeled.spp$Port.Match[1]<- paste(landings.prop.modeled.spp$Port.Match[1], "NY", sep = "")
landings.prop.modeled.spp$Port.Match[29]<- paste(landings.prop.modeled.spp$Port.Match[29], "NY", sep = "")
landings.prop.modeled.spp$Port.Match[34]<- paste(landings.prop.modeled.spp$Port.Match[34], "RI", sep = "")
landings.prop.modeled.spp$Port.Match[81]<- paste(landings.prop.modeled.spp$Port.Match[81], "RI", sep = "")

# Want to add species changes here, too....
port.sdm.dat<- read_csv(paste(out.dir, "PortData10142018.csv", sep = ""))

# Only need the all gear type...
port.sdm.dat.sub<- port.sdm.dat[port.sdm.dat$Gear == "All" & port.sdm.dat$Footprint == "Regular" & port.sdm.dat$Stat == "Mean", ]
port.sdm.dat.sub$Port.Match<- paste(port.sdm.dat.sub$Community_Only, port.sdm.dat.sub$State_Only, sep = ".")
port.sdm.dat.sub<- port.sdm.dat.sub %>%
  dplyr::select(., COMNAME, Port.Match, Baseline_sdm_p, Baseline_combo_b, Future_mean_sdm_p, Future_cold_sdm_p, Future_warm_sdm_p, Future_mean_combo_b, Future_cold_combo_b, Future_warm_combo_b, Future_mean_diff_sdm_p, Future_mean_diff_combo_b, Future_cold_diff_sdm_p, Future_cold_diff_combo_b, Future_warm_diff_sdm_p)

# Calculate percent changes...
perc_change<- function(base, future, set.zero = 0){
  base.adj<- ifelse(base == 0, set.zero, base)
  perc.diff<- 100*((future - base.adj)/base.adj)
  if(is.infinite(perc.diff)){
    perc.diff<- 5000
  } else {
    perc.diff
  }
  return(perc.diff)
}

port.sdm.dat.sub<- port.sdm.dat.sub %>%
  mutate(., "Future_mean_percdiff_sdm_p" = as.numeric(map2(Baseline_sdm_p, Future_mean_sdm_p, possibly(perc_change, NA))),
         "Future_cold_percdiff_sdm_p" = as.numeric(map2(Baseline_sdm_p, Future_cold_sdm_p, possibly(perc_change, NA))),
         "Future_warm_percdiff_sdm_p" = as.numeric(map2(Baseline_sdm_p, Future_warm_sdm_p, possibly(perc_change, NA))),
         "Furure_mean_percdiff_combo_b" = as.numeric(map2(Baseline_combo_b, Future_mean_combo_b, possibly(perc_change, NA))),
         "Furure_cold_percdiff_combo_b" = as.numeric(map2(Baseline_combo_b, Future_cold_combo_b, possibly(perc_change, NA))),
         "Furure_warm_percdiff_combo_b" = as.numeric(map2(Baseline_combo_b, Future_warm_combo_b, possibly(perc_change, NA))))

# Join...
landings.sdmchanges.port.out<- landings.prop.modeled.spp %>%
  dplyr::select(., COMNAME, Port.long, Port.Match, Landings, Prop.Landings, DevExp.P, DevExp.B, AUC.SDM, Calib.SDM, RMSE.SDM.P, RMSE.NEVA.B) %>%
  left_join(., port.sdm.dat.sub, by = c("COMNAME", "Port.Match"))
write.csv(landings.sdmchanges.port.out, file = paste(out.dir, "Species_Port_Landings_DistChangesforAllGears.csv"))

# Join to modeled species and sum prop.landings to get aggregate across species at a port level
landings.prop.modeled<- landings.prop %>%
  left_join(., dat.mod.results, by = "Match.Col") %>%
  drop_na(., AUC.SDM) %>%
  group_by(., Port.long) %>%
  summarize_at(., "Prop.Landings", sum, na.rm = T)
landings.prop.modeled$Port.Match<- as.character(landings.prop.modeled$Port.long)
landings.prop.modeled$Port.Match[1]<- paste(landings.prop.modeled$Port.Match[1], "NY", sep = "")
landings.prop.modeled$Port.Match[29]<- paste(landings.prop.modeled$Port.Match[29], "NY", sep = "")
landings.prop.modeled$Port.Match[34]<- paste(landings.prop.modeled$Port.Match[34], "RI", sep = "")
landings.prop.modeled$Port.Match[81]<- paste(landings.prop.modeled$Port.Match[81], "RI", sep = "")

# Want to add species changes here, too....
port.sdm.dat<- read_csv(paste(out.dir, "PortData10142018.csv", sep = ""))

# Only need the all gear type...
port.sdm.dat.sub<- port.sdm.dat[port.sdm.dat$Gear == "All" & port.sdm.dat$Footprint == "Regular" & port.sdm.dat$Stat == "Mean", ]
port.sdm.dat.sub$Landings.Match<- paste(port.sdm.dat.sub$Community_Only, port.sdm.dat.sub$State_Only, sep = ".")


# Plot?
ports.geo<- read_csv("./Data/geocodedports.csv")
ports.geo$Port.Match<- paste(gsub(" ", "_", ports.geo$PORT_NAME), ".", ports.geo$STATE, sep = "")
landings.prop.modeled<- landings.prop.modeled %>%
  left_join(., ports.geo, by = "Port.Match")

## Libraries you may need to install these...
library(tidyverse)
library(raster)
library(rgeos)
library(ggplot2)
library(viridis)

## Load in your dataset, at a minimum needs description, lat, long, value
dat<- landings.prop.modeled

# Spatial stuff -- gets us the states and shoreline
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
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

# Alright, plot time
dat<- dat %>%
  data.frame()


plot.out<- ggplot() + 
  # Update "fill" and "color" to change map color
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "#d9d9d9", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "#d9d9d9", color = "gray45", size = 0.15) +
  geom_point(data = dat, aes(x = lon, y = lat, fill = dat[,which(colnames(dat) == "Prop.Landings", arr.ind = T)]), shape = 21, size = 4.5, alpha = 0.65) +
  # Here you'd make adjustments to the point colors...
  scale_fill_gradient2(name = "Modeled species\n proportion of landings", low = "blue", mid = "white", high = "red", midpoint = 0.5, limits = c(0,1)) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"))

temp<- landings.prop.modeled %>%
  arrange(., -Prop.Landings, Port.long)
temp


# Only ports with 0.75 or greater prop landings
dat.sub<- dat %>%
  filter(., Prop.Landings >= 0.75)

plot.out<- ggplot() + 
  # Update "fill" and "color" to change map color
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "#d9d9d9", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "#d9d9d9", color = "gray45", size = 0.15) +
  geom_point(data = dat.sub, aes(x = lon, y = lat), shape = 21, size = 3, fill = "#31a354") +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"))

ggsave("~/Desktop/PortsGreaterThan0.75PropLandings.jpg", plot.out, height = 11, width = 8)

