quant_qual_plots_func<- function(quant.dat, qual.dat) {
  
  # Details -----------------------------------------------------------
  # The function plots our quantitative projections against Jon Hare's classification
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "gridExtra")
  package_check(packages.needed)
  
  if(FALSE) {
    quant.dat = "./Data/sdm.projections.rds"
    qual.dat = "./Data/neva.boot.results.csv"
  }
  
  # Read in the datasets
  sdm.proj<- readRDS(quant.dat)
  neva.dat<- read_csv(qual.dat)
  
  # First plot: Shelfwide average difference by NEVA directional effect, with color corresponding to the NEVA vulnerability rank
  # Getting shelfwide average
  # Function to get average shelfwide difference
  shelf_avg_difference<- function(df, scenario) {
    names(df)[4:6]<- c("X2025", "X2040", "X2055")
    diff<- switch(scenario,
           "2025" = df$X2025 - df$Baseline,
           "2040" = df$X2040 - df$Baseline,
           "2055" = df$X2055 - df$Baseline)
    return(round(mean(diff, na.rm = T), 3))
  }
  
  sdm.proj<- sdm.proj %>%
    mutate(., "diff.2025" = purrr::map2(Projections, "2025", possibly(shelf_avg_difference, NA)),
           "diff.2040" = purrr::map2(Projections, "2040", possibly(shelf_avg_difference, NA)),
           "diff.2055" = purrr::map2(Projections, "2055", possibly(shelf_avg_difference, NA)))
  
  # Bring in NEVA data
  neva<- read_csv(qual.dat) %>%
    dplyr::select(., COMNAME, HareRankVulnerability, HareDir.Eff)
  neva$COMNAME<- toupper(neva$COMNAME)
  
  # Join the SDM with NEVA for plottinf
  plot.one.dat.w<- sdm.proj %>%
    dplyr::select(., COMNAME, SEASON, Component, diff.2025, diff.2040, diff.2055) %>%
    unnest() %>%
    left_join(., neva, by = "COMNAME") %>%
    data.frame()
  
  # Move from wide to long format
  plot.one.dat<- plot.one.dat.w %>%
    gather(., Season.Time, Avg.Difference, -COMNAME, - SEASON, -Component, -HareRankVulnerability, -HareDir.Eff) %>%
    na.omit(.) 
  
  plot.one.dat$HareDir.Eff<- factor(plot.one.dat$HareDir.Eff, levels = c("Negative", "Neutral", "Positive"))
  plot.one.dat$HareRankVulnerability<- factor(plot.one.dat$HareRankVulnerability, levels = c("Low", "Moderate", "High", "Very High"))
  
  # Make the plots, season-component scatter plot
  plot.one.fall.p<- ggplot() +
    geom_point(data = subset(plot.one.dat, SEASON == "FALL" & Component == "Presence"), aes(x = HareDir.Eff, y = Avg.Difference, color = HareRankVulnerability), shape = 19, size = 2) +
    scale_color_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    facet_wrap(~Season.Time) +
    theme_dark() +
    theme(text = element_text(size = 16)) +
    ggtitle(paste("Fall.Presence"))
  
  plot.one.fall.b<- ggplot() +
    geom_point(data = subset(plot.one.dat, SEASON == "FALL" & Component == "Biomass"), aes(x = HareDir.Eff, y = Avg.Difference, color = HareRankVulnerability), shape = 19, size = 2) +
    scale_color_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    facet_wrap(~Season.Time) +
    theme_dark() +
    theme(text = element_text(size = 16)) +
    ggtitle(paste("Fall.Biomass"))
  
  plot.one.spring.p<- ggplot() +
    geom_point(data = subset(plot.one.dat, SEASON == "SPRING" & Component == "Presence"), aes(x = HareDir.Eff, y = Avg.Difference, color = HareRankVulnerability), shape = 19, size = 2) +
    scale_color_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    facet_wrap(~Season.Time) +
    theme_dark() +
    theme(text = element_text(size = 16)) +
    ggtitle(paste("Spring.Presence"))
  
  plot.one.spring.b<- ggplot() +
    geom_point(data = subset(plot.one.dat, SEASON == "SPRING" & Component == "Biomass"), aes(x = HareDir.Eff, y = Avg.Difference, color = HareRankVulnerability), shape = 19, size = 2) +
    scale_color_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    facet_wrap(~Season.Time) +
    theme_dark() +
    theme(text = element_text(size = 16)) +
    ggtitle(paste("Spring.Biomass"))
  
  # Panel plot layout
  png(file = "./ExploratoryResults/SDMAvgDiff.vs.NEVADirEff.png", width = 25, height = 8, units = "in", res = 100)
  grid.arrange(plot.one.fall.p, plot.one.fall.b, plot.one.spring.p, plot.one.spring.b, ncol = 2)
  dev.off()
  
  # Clearly some outlying species, which ones are these?
  # fall.p.df
  fall.p.df<- plot.one.dat %>%
    filter(SEASON == "FALL" & Component == "Presence", Season.Time == "diff.2055") %>%
    arrange(., HareRankVulnerability, HareDir.Eff, Avg.Difference)
  
  fall.p.spplot<- ggplot(data = fall.p.df, aes(x = HareDir.Eff, y = Avg.Difference, color = HareRankVulnerability)) + 
    scale_color_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    geom_text(aes(label=tolower(COMNAME)),size=3) +
    theme_dark()
  
  # spring.p.df
  spring.p.df<- plot.one.dat %>%
    filter(SEASON == "SPRING" & Component == "Presence", Season.Time == "diff.2055") %>%
    arrange(., HareRankVulnerability, HareDir.Eff, Avg.Difference)
  
  spring.p.spplot<- ggplot(data = spring.p.df, aes(x = HareDir.Eff, y = Avg.Difference, color = HareRankVulnerability)) + 
    scale_color_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    geom_text(aes(label=tolower(COMNAME)),size=3) +
    theme_dark()
  
  # Panel plot layout
  png(file = "./ExploratoryResults/SpeciesNamePlots.png", width = 12, height = 10, units = "in", res = 100)
  grid.arrange(fall.p.spplot, spring.p.spplot, nrow = 2)
  dev.off()
  
  # Box and whisker plot of p differences (adjustment) by directional effect, with bins for each vulnerability
  bw.fall<- ggplot(data = subset(plot.one.dat, SEASON == "FALL" & Component == "Presence"), aes(x = HareDir.Eff, y = Avg.Difference, fill = HareRankVulnerability)) +
    geom_boxplot() +
    stat_summary(fun.y = mean, color = "blue", geom = "point", shape = 18, size = 3.5, position=position_dodge(width=0.75)) +
    facet_wrap(~ Season.Time, nrow = 1) +
    scale_fill_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    theme_bw() 
  
  bw.spring<- ggplot(data = subset(plot.one.dat, SEASON == "SPRING" & Component == "Presence"), aes(x = HareDir.Eff, y = Avg.Difference, fill = HareRankVulnerability)) +
    geom_boxplot() +
    stat_summary(fun.y = mean, color = "blue", geom = "point", shape = 18, size = 3.5, position=position_dodge(width=0.75)) +
    facet_wrap(~ Season.Time, nrow = 1) +
    scale_fill_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    theme_bw() 
  
  png(file = "./ExploratoryResults/SDMBWPresenceAdjustment.vs.NEVADirEff.png", width = 12, height = 10, units = "in", res = 100)
  grid.arrange(bw.fall, bw.spring, nrow = 2)
  dev.off()
  
  # Third plot mean adjustments
  # First step, is need to get an average average difference per vulnerability rank and per directional effect -- all data, I think is in plot.one.dat
  # Group and average 
  adjustment.l <- plot.one.dat %>% 
    dplyr::filter(., Component == "Presence") %>%
    dplyr::group_by(., HareDir.Eff, HareRankVulnerability, SEASON, Season.Time) %>%
    dplyr::summarise(., "Adjustment.Mean" = mean(Avg.Difference, na.rm = T)) 
    
  count.l<- plot.one.dat %>% 
    dplyr::filter(., Component == "Presence") %>%
    dplyr::group_by(., HareDir.Eff, HareRankVulnerability, SEASON, Season.Time) %>%
    dplyr::summarise(., "Count" = n())
  
  plot.two.dat<- adjustment.l %>%
    dplyr::left_join(., count.l, by = c("HareDir.Eff", "HareRankVulnerability", "SEASON", "Season.Time"))
  
  fall.means<- ggplot(data = subset(plot.two.dat, SEASON == "FALL"), aes(x = HareDir.Eff, y = Adjustment.Mean, fill = HareRankVulnerability)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    geom_text(data = subset(plot.two.dat, SEASON == "FALL"), aes(label = Count, x = HareDir.Eff, vjust = ifelse(Adjustment.Mean >= 0, 0, 1)), color = "white", position = position_dodge(width = 0.8)) +
    facet_wrap(~ Season.Time, nrow = 1) +
    scale_fill_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    theme_dark() +
    ggtitle("Fall Means")
  
  spring.means<- ggplot(data = subset(plot.two.dat, SEASON == "SPRING"), aes(x = HareDir.Eff, y = Adjustment.Mean, fill = HareRankVulnerability)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    geom_text(data = subset(plot.two.dat, SEASON == "SPRING"), aes(label = Count, x = HareDir.Eff, vjust = ifelse(Adjustment.Mean >= 0, 0, 1)), color = "white", position = position_dodge(width = 0.8)) +
    facet_wrap(~ Season.Time, nrow = 1) +
    scale_fill_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
    theme_dark() +
    ggtitle("Spring Means")
  
  png(file = "./ExploratoryResults/SDMMeanAdjustment.vs.NEVADirEff.png", width = 14, height = 8, units = "in", res = 100)
  grid.arrange(fall.means, spring.means, nrow = 2)
  dev.off()
 
  # Third plot: NEVA Certainty in Vulnerability vs. SDM AUC
  # First need to get all the AUC data, one file per group
  # Loop over mod.table
  # files.list<- list.files(out.dir, ".mod.table.csv")
  # 
  # for(i in 1:length(files.list)) {
  #   file.use<- files.list[i]
  #   
  #   temp<- read.csv(paste(out.dir, file.use, sep = ""))
  #   temp.nona<- temp[!is.na(temp$Dev.Exp.p),]
  #   
  #   if(i == 1) {
  #     mod.tab.out<- temp.nona
  #   } else {
  #     mod.tab.out<- dplyr::bind_rows(mod.tab.out, temp.nona)
  #   }
  # }
  # 
  # # Check for duplicates...
  # spp.dups<- mod.tab.out$COMNAME[duplicated(mod.tab.out$COMNAME)]
  # dups<- mod.tab.out[mod.tab.out$COMNAME == spp.dups,]
  # mod.tab.out<- mod.tab.out[-as.numeric(row.names(dups[!complete.cases(dups),])),]
  # 
  # # Now need to merge this with the neva data
  # plot.three.dat<- neva %>%
  #   dplyr::select(., COMNAME, Vulnerability.Certainty, Dir.Eff.Certainty) %>%
  #   dplyr::left_join(., mod.tab.out, by = "COMNAME") %>%
  #   gather(., "Stat", "Value", -COMNAME, -Vulnerability.Certainty, -Dir.Eff.Certainty, -X)
  # 
  # plot.three.dat$Vulnerability.Certainty<- factor(plot.three.dat$Vulnerability.Certainty, levels = c("Low", "Moderate", "High", "Very High"))
  # plot.three.dat$Dir.Eff.Certainty<- factor(plot.three.dat$Dir.Eff.Certainty, levels = c("Low", "Moderate", "High", "Very High"))
  # 
  # plot.three.a<- ggplot() +
  #   geom_point(data = plot.three.dat, aes(x = Value, y = Vulnerability.Certainty)) +
  #   facet_wrap(~Stat, scales = "free")
  # 
  # plot.three.b<- ggplot() +
  #   geom_point(data = plot.three.dat, aes(x = Value, y = Dir.Eff.Certainty)) +
  #   facet_wrap(~Stat, scales = "free")
  # 
  # png(file = paste(out.dir, "NEVACertainty.vs.SDMAUC.png", sep = ""), width = 14, height = 8, units = "in", res = 100)
  # grid.arrange(plot.three.a, plot.three.b, ncol = 2)
  # dev.off()
  
  # Save a few of the outputs
  out.1<- data.frame(plot.one.dat)
  write.csv(out.1, file = "./Data/Species.DirEff.VulnRank.Diffs.csv")
  
  out.2<- data.frame(plot.two.dat)
  write.csv(out.2, file = "./Data/AverageAdjustment.DirEff.VulnRank.csv")
}