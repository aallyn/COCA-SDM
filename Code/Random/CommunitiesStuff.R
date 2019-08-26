##### Communities!!!!
###########

## Read in Footprint data 
all.foot.dat<- readRDS("~/GitHub/COCA/Data/VTR fishing footprints by community and gear type 2011-2015.rds")
ports.names<- all.foot.dat$JGS.COMMUNITY

safe.foots<- raster::stack(unlist(lapply(all.foot.dat$JGS.NOAA.SAFE.COMMUNITY.GEAR.FOOTPRINTS, "[", 3)))
print(safe.foots)

## Need better description of Community Names and Gear Types..
# Gear Types are in COST ID list item
gear.types<- all.foot.dat$COST_ID
gear.types<- ifelse(gear.types == 1, "Dredge",
                    ifelse(gear.types == 2, "Gillnet",
                           ifelse(gear.types == 3, "Longline",
                                  ifelse(gear.types == 4, "Pot/Trap",
                                         ifelse(gear.types == 5, "Purse/Seine",
                                                ifelse(gear.types == 6, "Trawl", "Other"))))))

# Community names in JGS.Community list item
safe.foot.names<- paste(all.foot.dat$JGS.COMMUNITY, gear.types, sep = "-")
names(safe.foots)<- safe.foot.names

# Plot them!
if(FALSE){
  
  ## Object renaming
  port.foots.stack<- safe.foots
  
  ## Ports to Plot
  focal.ports<- c("WESTPORT_MA", "CUSHING_ME")
  
  ## Some functions for messing around with the names
  gsub_func1<- function(pattern, replacement, string){
    out<- gsub(paste(pattern, ".", sep = ""), replacement, string)
    return(out)
  }
  
  gsub_func2<- function(pattern, replacement, string){
    out<- gsub(pattern, replacement, string)
    return(out)
  }
  
  ## Spatial stuff
  #Bounds
  xlim.use<- c(-77, -65)
  ylim.use<- c(35, 45)
  
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
  
  port.foots.focal<-  port.foots.stack[[which(str_detect(names(port.foots.stack), paste(focal.ports, collapse = "|")))]]
  gear.types<- c("Dredge", "Gillnet", "Longline", "Pot.Trap", "Purse.Seine", "Trawl", "Other")
  
  # Make panel plot for each gear type: Ports as rows, Fishing footprint as columns
  for(i in seq_along(gear.types)){
    if(i == 5){
      next()
    }
    port.foots.use<- port.foots.focal[[which(str_detect(names(port.foots.focal), gear.types[i]))]]
    names(port.foots.use)<- str_replace_all(names(port.foots.use), c(".JGS.SAFE.PROPORTION" = ""))
    
    out.df<- raster::as.data.frame(port.foots.use, xy = T) %>%
      gather(., "Port_Gear", "Z", -x, -y)
    
    gear.replace<- paste('.', gear.types[i], sep = '')
    out.df$Port<- gsub(gear.replace, "", out.df$Port_Gear)
    out.df$Port<- gsub("MaxD", "", out.df$Port)
    out.df$Port<- gsub("1.5x", "", out.df$Port)
    out.df<- out.df %>%
      as_tibble() %>%
      mutate(., "Gear" = pmap(list(pattern = Port, replacement = list(""), string = Port_Gear), gsub_func1)) %>%
      data.frame()
    out.df$Footprint<- ifelse(str_detect(out.df$Port_Gear, "1.5x"), "1.5xExtended", 
                              ifelse(str_detect(out.df$Port_Gear, "MaxD"), "Extended", "Regular"))
    out.df$Z<- ifelse(out.df$Z > 0, "Fished", out.df$Z)
    out.df$Footprint<- factor(out.df$Footprint, levels = c("Regular", "Extended", "1.5xExtended"))
    
    plot.out<- ggplot() + 
      geom_tile(data = out.df, aes(x = x, y = y, fill = Z), show.legend = FALSE) +
      scale_fill_manual(values = "#31a354", na.value = "white") +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.text=element_text(size=8), legend.title=element_text(size=8)) +
      ggtitle(gear.types[i]) +
      facet_wrap(~Port+Footprint, ncol = 3)
    ggsave(paste("~/GitHub/COCA/Data/", gear.types[i], "Footprints.jpg", sep = ""), plot.out, width = 11, height = 8, units = "in")
    print(i)
  }
}

length(port.names)
