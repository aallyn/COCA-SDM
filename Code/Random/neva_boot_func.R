neva_boot_func<- function(qual.sensitivity.exposure.path = "./Data/JHareQualitativeDataResults.csv", qual.directional.effect.path = "./Data/JHareDirectionalEffect.csv", boot.run = 10000) {
  
  # Details -----------------------------------------------------------
  # The function reads in Hare et al voting dataset. The original dataset is in a wide format, so the function then transforms it to a long/tidy version, such that each row is a species-attribute-vote. The function then goes through the process implemented by Hare et al. to calculate an overall vulnerability rank, directional effect and uncertainty in both of those. 
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse")
  package_check(packages.needed)
  
  if(FALSE) {
    qual.sensitivity.exposure.path = "./Data/JHareQualitativeDataResults.csv"
    qual.directional.effect.path = "./Data/JHareDirectionalEffect.csv"
    boot.run = 1000
    prob.spp = c("American Lobster", "Winter Skate", "Spiny Dogfish")
  }
  
  # Need to create species adjustment dataset by bootstrapping voting dataset
  # First need to rearrange qualitative data
  # Sensitivity/Exposure
  qual.dat<- read_csv(qual.sensitivity.exposure.path)
  
  # Wide to long format
  qual.dat.l<- gather(qual.dat, "Score", "Votes", Low, Moderate, High, Very.High)
  qual.dat.l<- arrange(qual.dat.l, Species, Functional.Group, Attribute, Score)
  
  # For bootstrapping, really need each row as a vote....
  qual.dat.boot<- qual.dat.l[rep(seq.int(1,nrow(qual.dat.l)), qual.dat.l$Votes), 1:5]
  qual.dat.boot$id<- seq(1, nrow(qual.dat.boot), 1)
  
  # Directional effect
  dir.eff.dat<- read.csv(qual.directional.effect.path)
  dir.eff.dat.l<- gather(dir.eff.dat, "Score", "Votes", Negative, Neutral, Positive)
  dir.eff.dat.l<- arrange(dir.eff.dat.l, Species, Functional.Group, Score)
  
  # For bootstrapping, really need each row as a vote....
  dir.eff.dat.boot<- dir.eff.dat.l[rep(seq.int(1,nrow(dir.eff.dat.l)), dir.eff.dat.l$Votes), 1:4]
  dir.eff.dat.boot$id<- seq(1, nrow(dir.eff.dat.boot), 1)
  
  # Bootstrap result allocation
  # Qualitative data
  qual.boot.res<- data.frame(matrix(nrow = length(unique(qual.dat.boot$Species)), ncol = boot.run+1))
  names(qual.boot.res)<- c("COMNAME")
  qual.boot.res$COMNAME<- unique(qual.dat.boot$Species)
  
  # Directional effect
  dir.eff.boot.res<- data.frame(matrix(nrow = length(unique(dir.eff.dat.boot$Species)), ncol = boot.run+1))
  names(dir.eff.boot.res)<- c("COMNAME")
  dir.eff.boot.res$COMNAME<- unique(dir.eff.dat.boot$Species)
  
  for(i in 1:boot.run) {
    
    ## Sensitivity and Exposure Data
    # Get bootstrapped sample
    dat.temp<- qual.dat.boot %>%
      group_by(., Species, Attribute) %>%
      do(sample_n(., nrow(.), replace = T)) %>%
      arrange(., Species, Functional.Group, Attribute, Score) %>%
      data.frame(.)
    
    # Aggregate to calculate weighted averages
    # Aggregate
    dat.agg.temp<- dat.temp %>%
      group_by(., Species, Functional.Group, Attribute, Attribute.Category, Score) %>%
      dplyr::summarise(Votes = n()) %>%
      data.frame(.)
    
    # Convert Score to numeric (1:4)
    dat.agg.temp$Score<- factor(dat.agg.temp$Score, levels = c("Low", "Moderate", "High", "Very.High"))
    dat.agg.temp$Score.Numeric<- as.numeric(dat.agg.temp$Score)
    
    # Calculate weighted mean by species-attribute-attribute.category, then characterize these into low, moderate, high and very high. This will allow us to count them and apply the logic rule used by Hare et al.
    dat.wt.mean<- dat.agg.temp %>%
      group_by(., Species, Attribute, Attribute.Category) %>%
      dplyr::summarise("Weighted.Mean" = weighted.mean(Score.Numeric, Votes)) %>%
      dplyr::mutate(., "Logic.Rule" = cut(Weighted.Mean, breaks = c(0, 2.5, 3.0, 3.5, 100), labels = c("Low", "Moderate", "High", "Very.High")))
    
    # Get the counts (number of votes in each species - attribute category - vulneraibity) to apply logic rule              
    dat.scores<- dat.wt.mean %>%
      group_by(., Species, Attribute.Category, Logic.Rule) %>%
      dplyr::summarise("Number" = n()) %>%
      data.frame(.)
    
    # Apply logic rule to translate number of votes to a low-moderatre-high-very high rank for each species - attribute category
    spp.ranks.split<- dat.scores %>%
      dplyr::mutate(Rank.Temp = as.numeric(ifelse(Logic.Rule == "Very.High" & Number >= 3, 4, 
                                                  ifelse(Logic.Rule == "High" | Logic.Rule == "Very.High" & Number >= 2, 3, 
                                                         ifelse(Logic.Rule == "High" | Logic.Rule == "Very.High" | Logic.Rule == "Moderate" & Number >= 2, 2, 1))))) %>%
      dplyr::group_by(., Species, Attribute.Category) %>%
      dplyr::summarise(Rank = max(Rank.Temp)) 
    
    # Calculate one overall vulnerability rank for each species, which is exposure factor rank * sensitivity attribute rank (ranges 1-16)
    spp.ranks.overall<- spp.ranks.split %>%
      spread(., Attribute.Category, Rank) %>% 
      dplyr::group_by(., Species) %>%
      dplyr::summarise(Overall.Rank = Exposure.Factor*Sensitivity.Attribute) %>% 
      dplyr::mutate(., "Overall.Rank.Code" = cut(Overall.Rank, breaks = c(0, 3, 6, 9, 16), labels = c("Low", "Moderate", "High", "Very.High")))
    
    qual.boot.res[,i+1]<- spp.ranks.overall$Overall.Rank.Code
    colnames(qual.boot.res)[i+1]<- paste("boot.", i, sep = "")
    
    ## Directional Effect Data
    # Get bootstrapped sample
    dat.temp<- dir.eff.dat.boot %>%
      group_by(., Species) %>%
      do(sample_n(., nrow(.), replace = T)) %>%
      arrange(., Species, Score) %>%
      data.frame(.)
    
    # Aggregate to calculate weighted averages
    # Aggregate
    dat.agg.temp<- dat.temp %>%
      group_by(., Species, Score) %>%
      dplyr::summarise(Votes = n()) %>%
      data.frame(.)
    
    # Convert directional effect to number
    dat.agg.temp$Score.Numeric<- ifelse(dat.agg.temp$Score == "Negative", -1,
                                        ifelse(dat.agg.temp$Score == "Neutral", 0, 1))
    
    
    # Calculate weighted mean by species, then characterize these into negative. neutral, positive. 
    dat.wt.mean<- dat.agg.temp %>%
      group_by(., Species) %>%
      dplyr::summarise("Weighted.Mean" = weighted.mean(Score.Numeric, Votes)) %>%
      dplyr::mutate(., "Directional.Effect" = cut(Weighted.Mean, breaks = c(-100, -0.333, 0.333, 100), labels = c("Negative", "Neutral", "Positive")))
    
    dir.eff.boot.res[,i+1]<- dat.wt.mean$Directional.Effect
    colnames(dir.eff.boot.res)[i+1]<- paste("boot.", i, sep = "")
  }
  
  # Now have the overall vulnerability rank and the directional effect for each bootstrap. According to Hare et al. (2016), we then want to calculate the proportion of boot runs in each of the vulnerability rank bins and in each of the directional effect bins, by species.
  # Vulnerability first
  # Currently, data is in the wide format, I think for this it would be easier to put things in the long format and then calculate the summary per species-vulnearbility rank or per species-directional effect.
  qual.boot.res.l<- qual.boot.res %>%
    tidyr::gather(., Boot.Run, Rank, -COMNAME)
  
  # Proportions
  qual.boot.res.props<- qual.boot.res.l %>%
    dplyr::group_by(COMNAME, Rank) %>% 
    summarise("Count" = n()) %>%
    mutate(., "Prop" = Count/boot.run)
  
  # Final dataset -- keep maximum proportion, and then assign certainty based on Hare et al. (2016), with very high certainty (prop >95%), high certainty (prop 90-95%), moderate certainty (66-90%), low certainty (<66%).
  qual.res.final<- qual.boot.res.props %>%
    dplyr::group_by(COMNAME) %>%
    dplyr::filter(Prop == max(Prop)) %>%
    dplyr::mutate(., "Vulnerability.Certainty" = cut(Prop, breaks = c(0, 0.66, 0.90, 0.95, 1.0), labels = c("Low", "Moderate", "High", "Very High")))
  
  # Directional effect next, same steps
  # Wide to long formate
  dir.eff.boot.res.l<- dir.eff.boot.res %>%
    tidyr::gather(., Boot.Run, Rank, -COMNAME)
  
  # Proportions
  dir.eff.boot.res.props<- dir.eff.boot.res.l %>%
    dplyr::group_by(COMNAME, Rank) %>% 
    summarise("Count" = n()) %>%
    mutate(., "Prop" = Count/boot.run)
  
  # Final dataset
  dir.eff.res.final<- dir.eff.boot.res.props %>%
    dplyr::group_by(COMNAME) %>%
    dplyr::filter(Prop == max(Prop)) %>%
    dplyr::mutate(., "Dir.Eff.Certainty" = cut(Prop, breaks = c(0, 0.66, 0.90, 0.95, 1.0), labels = c("Low", "Moderate", "High", "Very High")))
  
  ## Join and write out the results
  neva.out<- qual.res.final %>%
    dplyr::left_join(., dir.eff.res.final, by = "COMNAME") %>%
    data.frame()
  names(neva.out)<- c("COMNAME", "Vulnerability.Rank", "Vulnerability.Count", "Vulnerability.Prop", "Vulnerability.Certainty", "Dir.Eff", "Dir.Eff.Count", "Dir.Eff.Prop", "Dir.Eff.Certainty")
  write.csv(neva.out, file = paste(out.dir, "neva.boot.results.csv", sep = ""))
}