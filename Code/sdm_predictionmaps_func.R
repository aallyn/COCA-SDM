plot_func<- function(response, mod.fitted.p, mod.fitted.b, season, base.preds, fut.preds, like.posts, posts.samps) {
  if(FALSE){
    response = "Presence"
    mod.fitted.p = dat.use$mod.fitted.p[[1]]
    mod.fitted.b = dat.use$mod.fitted.b[[1]]
    season = season.use
    base.preds = base.preds
    fut.preds = fut.preds
    like.posts = like.posts.b
    posts.samps = posts.df.b
  }
  
  if(response == "Presence"){
    # SDM
    base.data<- base.preds$Data[[match(season, base.preds$SEASON)]]
    base.mu<- base.data %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
    base.map<- data.frame("x" = base.mu$x, "y" = base.mu$y, "pred" = predict.gam(mod.fitted.p, newdata = base.mu, type = "response"))
    
    # Mean climate model SST
    newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
    newdat.mu<- newdat.mu %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
    map.fut.mu<- data.frame("x" = newdat.mu$x, "y" = newdat.mu$y, "pred" = predict.gam(mod.fitted.p, newdata = newdat.mu, type = "response"))
    
    # pct05 climate model
    newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
    newdat.05<- newdat.05 %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
    names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
    map.fut.pct05<- data.frame("x" = newdat.05$x, "y" = newdat.05$y, "pred.05" = predict.gam(mod.fitted.p, newdata = newdat.05, type = "response"))
    
    # pct95 climate model
    newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
    newdat.95<- newdat.95 %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
    names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
    map.fut.pct95<- data.frame("x" = newdat.95$x, "y" = newdat.95$y, "pred.95" = predict.gam(mod.fitted.p, newdata = newdat.95, type = "response"))
    
    # Maps
    # Spatial projections
    proj.wgs84<- "+init=epsg:4326" #WGS84
    proj.utm<- "+init=epsg:2960" #UTM 19
    
    # NELME domaine
    nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
    st_crs(nelme)<- proj.wgs84
    nelme.sp<- as(nelme, "Spatial")
    
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
    
    # Ideally, want a bit more fine scale of a map for visual purposes. To do that, we will be doing some interpolation and then cropping the map to the NELME border. This can take a bit of time, so instead of doing it for EVERY dataset, lets do it once here so we can get a vector of the points to keep. Then just pass that into the plotting function.
    coords.df<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y)
    pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Baseline
    data.use<- base.map
    pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)

    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")

    # Discrete scale
    pred.df.use$breaks<- cut(pred.df.use$z,
                       breaks = seq(from = 0.0, to = 1.0, by = 0.1),
                       labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
    pred.df.use<- pred.df.use %>%
      drop_na(., breaks)

    plot.base.sdm<- ggplot() +
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "SDM Base", na.value = "white", discrete = T, drop = FALSE, direction = 1) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) +
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) +
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) +
      guides(fill=guide_legend(ncol=2))
    
    # Future mean
    data.use<- map.fut.mu
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    # Discrete scale
    pred.df.use$breaks<- cut(pred.df.use$z, 
                             breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                             labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
    pred.df.use<- pred.df.use %>%
      drop_na(., breaks)
    
    plot.fut<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "SDM P(Pres)\nAverage future", na.value = "white", discrete = T, drop = F) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) +
      guides(fill=guide_legend(ncol=2))
    
    # Lower
    data.use<- map.fut.pct05
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred.05))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    # Discrete scale
    pred.df.use$breaks<- cut(pred.df.use$z, 
                             breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                             labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
    pred.df.use<- pred.df.use %>%
      drop_na(., breaks)
    
    plot.fut.lwr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "SDM P(Pres)\nCold future", na.value = "white", discrete = T, drop = F) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) +
      guides(fill=guide_legend(ncol=2))
    
    # Upper
    data.use<- map.fut.pct95
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred.95))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    # Discrete scale
    pred.df.use$breaks<- cut(pred.df.use$z, 
                             breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                             labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
    pred.df.use<- pred.df.use %>%
      drop_na(., breaks)
    
    plot.fut.upr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "SDM P(Pres)\nWarm future", na.value = "white", discrete = T, drop = F) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) +
      guides(fill=guide_legend(ncol=2))
    
    # Differences
    # Getting limits
    # sdm<- range(data.frame("x" = fut.map$x, "y" = fut.map$y, "pred" = fut.map$pred - base.map$pred)$pred, na.rm = T)
    sdm<- data.frame("x" = map.fut.mu$x, "y" = map.fut.mu$y, "pred" = map.fut.mu$pred - base.map$pred)
    sdm.range<- range(sdm$pred, na.rm = T)
    sdm.lwr<- data.frame("x" = map.fut.pct05$x, "y" = map.fut.pct05$y, "pred" = map.fut.pct05$pred.05 - base.map$pred)
    sdm.lwr.range<- range(sdm.lwr$pred, na.rm = T)
    sdm.upr<- data.frame("x" = map.fut.pct95$x, "y" = map.fut.pct95$y, "pred" = map.fut.pct95$pred.95 - base.map$pred)
    sdm.upr.range<- range(sdm.upr$pred, na.rm = T)
    diff.lim<- c(round(min(sdm.range, sdm.lwr.range, sdm.upr.range), 2), round(max(sdm.range, sdm.lwr.range, sdm.upr.range), 2))
    
    # SDM
    data.use<- sdm
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.diff<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(name = "SDM P(Pres) Diff\nAverage future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 
    
    # Lower
    data.use<- sdm.lwr
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.diff.lwr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(name = "SDM P(Pres) Diff\nCold future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 
    
    # Lower
    data.use<- sdm.upr
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.diff.upr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(name = "SDM P(Pres) Diff\nWarm future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 
    
    out<- plot_grid(plot.base.sdm, NULL, NULL, plot.fut, plot.fut.lwr, plot.fut.upr, plot.diff, plot.diff.lwr, plot.diff.upr, nrow = 3, align = "hv", scale = 1.1)
  }
  
  if(response == "Biomass"){
    ilink <- family(mod.fitted.b)$linkinv
    # Get SDM predictions
    # Get maximimum value and extract row from posts.samps
    
    mod.ind.b<- like.posts.b$Iteration[which.max(like.posts.b$Value)]
    best.fit.b<- posts.samps.b[mod.ind.b,]
    best.fit.mat.b<- matrix(as.numeric(best.fit.b), nrow = 1, ncol = length(best.fit.b), byrow = T, dimnames = list(NULL, names(coef(dat.use$mod.fitted.b))))
    
    # SDM
    sdm.base<- base.preds$Data[[match(season, base.preds$SEASON)]]
    sdm.base<- sdm.base %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
    sdm.base<- predict.gam(mod.fitted.p, newdata = sdm.base, type = "response")
      
    
    # Mean climate model SST
    newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
    newdat.mu<- newdat.mu %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
    sdm.fut.mu<- predict.gam(mod.fitted.p, newdata = newdat.mu, type = "response")
 
    # pct05 climate model
    newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
    newdat.05<- newdat.05 %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
    names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
    sdm.fut.pct05<- predict.gam(mod.fitted.p, newdata = newdat.05, type = "response")
 
    # pct95 climate model
    newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
    newdat.95<- newdat.95 %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
    names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
    sdm.fut.pct95<- predict.gam(mod.fitted.p, newdata = newdat.95, type = "response")
   
    # Getting combined biomass predictions
    # Get maximimum value and extract row from posts.samps
    mod.ind<- like.posts$Iteration[which.max(like.posts$Value)]
    best.fit<- posts.samps[mod.ind,]
    best.fit.mat<- matrix(as.numeric(best.fit), nrow = 1, ncol = length(best.fit), byrow = T, dimnames = list(NULL, names(coef(mod.fitted.b))))
    
    # Make predictions with these values
    lpmat.base<- predict.gam(mod.fitted.b, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "lpmatrix")
    combo.map.base<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base * exp(ilink(as.numeric(lpmat.base %*% t(best.fit.mat)))))
    
    # Mean climate model SST
    newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
    newdat.mu<- newdat.mu %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
    lpmat.fut.mu<- predict.gam(mod.fitted.b, newdata = newdat.mu, type = "lpmatrix")
    combo.map.fut.mu<- data.frame("x" = newdat.mu$x, "y" = newdat.mu$y, "pred" = sdm.fut.mu * exp(ilink(as.numeric(lpmat.fut.mu %*% t(best.fit.mat)))))
    
    # pct05 climate model
    newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
    newdat.05<- newdat.05 %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.OISST.Scale"))
    names(newdat.05)[4]<- "SEASONALMU.OISST.Scale"
    lpmat.fut.pct05<- predict.gam(mod.fitted.b, newdata = newdat.05, type = "lpmatrix")
    combo.map.fut.pct05<- data.frame("x" = newdat.05$x, "y" = newdat.05$y, "pred.05" = sdm.fut.pct05 * exp(ilink(as.numeric(lpmat.fut.pct05 %*% t(best.fit.mat)))))
    
    # pct95 climate model
    newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
    newdat.95<- newdat.95 %>%
      unnest() %>%
      dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.OISST.Scale"))
    names(newdat.95)[4]<- "SEASONALMU.OISST.Scale"
    lpmat.fut.pct95<- predict.gam(mod.fitted.b, newdata = newdat.95, type = "lpmatrix")
    combo.map.fut.pct95<- data.frame("x" = newdat.95$x, "y" = newdat.95$y, "pred.95" = sdm.fut.pct95 * exp(ilink(as.numeric(lpmat.fut.pct95 %*% t(best.fit.mat)))))
    
    # Maps
    # Spatial projections
    proj.wgs84<- "+init=epsg:4326" #WGS84
    proj.utm<- "+init=epsg:2960" #UTM 19
    
    # NELME domaine
    nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
    st_crs(nelme)<- proj.wgs84
    nelme.sp<- as(nelme, "Spatial")
    
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
    
    # Ideally, want a bit more fine scale of a map for visual purposes. To do that, we will be doing some interpolation and then cropping the map to the NELME border. This can take a bit of time, so instead of doing it for EVERY dataset, lets do it once here so we can get a vector of the points to keep. Then just pass that into the plotting function.
    coords.df<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y)
    pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Combo
    data.use<- combo.map.base
    pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    # Limits
    combo.range<- range(combo.map.base$pred, na.rm = T)
    combo.fut.range<- range(combo.map.fut.mu$pred, na.rm = T)
    combo.lwr.range<- range(combo.map.fut.pct05$pred.05, na.rm = T)
    combo.upr.range<- range(combo.map.fut.pct95$pred.95, na.rm = T)
    combo.lim<- c(round(min(combo.range,  combo.fut.range, combo.lwr.range, combo.upr.range), 2), round(max(combo.range, combo.fut.range, combo.lwr.range, combo.upr.range), 2))
    
    # Continuous scale plot
    plot.base.combo<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nBase", na.value = "white", limits = combo.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))
    
    ## Future
    # Combo future
    data.use<- combo.map.fut.mu
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.fut.combo<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nAverage future", na.value = "white", limits = combo.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))
    
    # Lower
    data.use<- combo.map.fut.pct05
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")

    plot.fut.combo.lwr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nColder future", na.value = "white", limits = combo.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))
    
    # Upper
    data.use<- combo.map.fut.pct95
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")

    plot.fut.combo.upr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_viridis(option = "viridis", name = "NEVA Biomass\nWarmer future", na.value = "white", limits = combo.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      #geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      # geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8))
    
    # Differences
    # Getting limits
    combo<- data.frame("x" = combo.map.fut.mu$x, "y" = combo.map.fut.mu$y, "pred" = combo.map.fut.mu$pred - combo.map.base$pred)
    combo.range<- range(combo$pred, na.rm = T)
    combo.lwr<- data.frame("x" = combo.map.fut.pct05$x, "y" = combo.map.fut.pct05$y, "pred" = combo.map.fut.pct05$pred - combo.map.base$pred)
    combo.lwr.range<- range(combo.lwr$pred, na.rm = T)
    combo.upr<- data.frame("x" = combo.map.fut.pct95$x, "y" = combo.map.fut.pct95$y, "pred" = combo.map.fut.pct95$pred - combo.map.base$pred)
    combo.upr.range<- range(combo.upr$pred, na.rm = T)
    diff.lim<- c(round(min(combo.range, combo.lwr.range, combo.upr.range), 2), round(max(combo.range, combo.lwr.range, combo.upr.range), 2))
    
    # Combo
    data.use<- combo
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.diff.combo<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(name = "NEVA Biomass Diff\nAverage future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 
    
    # Lower
    data.use<- combo.lwr
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.diff.combo.lwr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(name = "NEVA Biomass Diff\nCold future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 
    
    # Lower
    data.use<- combo.upr
    pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    plot.diff.combo.upr<- ggplot() + 
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(name = "NEVA Biomass Diff\nWarm future", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("Lat") +
      scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=8), legend.title=element_text(size=8)) 
    
    out<- plot_grid(plot.base.combo, NULL, NULL, plot.fut.combo, plot.fut.combo.lwr, plot.fut.combo.upr, plot.diff.combo, plot.diff.combo.lwr, plot.diff.combo.upr, nrow = 3, align = "hv", scale = 1.1)
  }
  return(out)
}
