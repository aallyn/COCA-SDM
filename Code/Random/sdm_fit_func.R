########## A function for fitting SDMs to NOAA NEFSC cleaned boottom trawl survey data ##########
sdm_func<- function(model.dat = "./Data/model.dat.rds", model.framework = "GAM.DLN", split.type = "transfer", train.start = NULL, train.end = NULL, test.start = NULL, test.end = NULL) {
  
  # Details -----------------------------------------------------------
  # The function reads in the .rds data file, created by the trawl_data_prep function. After loading the file, the function then fits SDMs to each of the species in the NEVA assessment. For the fitting process, the user can select between a "GAM.DLN" model (generalized additive model delta lognormal), or a "RF.DLN" (random forest delta lognormal) model framework. In addition, the user can select how the dataset is split into training and testing pieces, by either selecting "transfer" and then specifying start and end dates for each component or by selecting "random" and then a 75% random subset of the data is used for training and the other 25% is used for testing. After fitting the models to each of the species, the function outputs a few different results. 
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "mgcv", "randomForest", "Boruta", "sp", "raster", "geosphere")
  package_check(packages.needed)
  
  # Debugging ---------------------------------------------------------
  if(FALSE) {
    model.dat<- "./Data/model.dat.rds"
    model.framework<- "GAM.DLN"
    split.type<- "transfer"
    train.start<- "1982-01-01"
    train.end<- "2010-12-31"
    test.start<- "2011-01-01"
    test.end<- "2016-01-01"
  }
  
  # Read in model data ---------------------------------------------------------
  # Read in the data, subset to only species in assessment
  fish.spp<- read.csv("./Data/Assesmentfishspecies.csv")
  dat<- readRDS(model.dat) %>% 
    filter(., SVSPP %in% fish.spp$SVSPP) %>%
    left_join(., fish.spp, by = "SVSPP")
  
  # Split data into testing and training datasets based on the splits ratio or dates
  if(split.type == "random") {
    n<- nrow(dat)
    groups<- rep(c("TRAIN", "TEST"), n*c(split, 1-split))
    dat$TRAIN.TEST<- sample(groups)
  } 
  
  if(split.type == "transfer") {
    dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                            ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither")) 
  }
  
  # Create nested dataframes, one for testing, one for training
  # Training
  dat.train<- dat %>%
    group_by(., COMNAME, SEASON, TRAIN.TEST) %>%
    dplyr::filter(., TRAIN.TEST == "TRAIN") %>%
    nest(.key = "TRAIN.DATA") %>%
    arrange(COMNAME)
  
  # Testing
  dat.test<- dat %>%
    group_by(., COMNAME, SEASON, TRAIN.TEST) %>%
    dplyr::filter(., TRAIN.TEST == "TEST") %>%
    nest(.key = "TEST.DATA") %>%
    arrange(COMNAME) 
  
  # Fit models -------------------------------------------------------
  # Model functions -- Create a seperate model function function with these in it?
  if(model.framework == "GAM.DLN") {
    m1.p <- function(df) {
      gam(PRESENCE ~ s(SHELF_POS, k = 4, fx = FALSE, bs = 'tp') + s(DEPTH, k = 4, fx = FALSE, bs = 'tp') + s(SEASONALMU.OISST, k = 4, fx = FALSE, bs = 'tp') + ti(DEPTH, SEASONALMU.OISST, k = 4, bs = c("tp", "tp")), data = df, family = binomial(link = logit), select = TRUE)
    }
    
    m1.b <- function(df) {
      gam(BIOMASS.MOD ~ s(SHELF_POS, k = 4, fx = FALSE, bs = 'tp') + s(DEPTH, k = 4, fx = FALSE, bs = 'tp') + s(SEASONALMU.OISST, k = 4, fx = FALSE, bs = 'tp') + ti(DEPTH, SEASONALMU.OISST, k = 4, bs = c("tp", "tp")), data = df, family = gaussian, select = TRUE)
    }
    
    # Temp deviance explained
    temp.dev.exp<- function(df, mod.fitted, response) {
      if(response == "PRESENCE"){
        null.mod<- gam(PRESENCE ~ 1, data = df)
        no.temp<- gam(PRESENCE ~ s(SHELF_POS, k = 4, fx = FALSE, bs = 'tp') + s(DEPTH, k = 4, fx = FALSE, bs = 'tp'), data = df, family = binomial(link = logit), sp = c(mod.fitted$sp[1], mod.fitted$sp[2], mod.fitted$sp[3], mod.fitted$sp[4]))
        temp.dev.exp<- (deviance(no.temp) - deviance(mod.fitted))/deviance(null.mod)
        return(temp.dev.exp)
      }
      
      if(response == "BIOMASS.MOD"){
        null.mod<- gam(BIOMASS.MOD ~ 1, data = df)
        no.temp<- gam(BIOMASS.MOD ~ s(SHELF_POS, k = 4, fx = FALSE, bs = 'tp') + s(DEPTH, k = 4, fx = FALSE, bs = 'tp'), data = df, family = gaussian, sp = c(mod.fitted$sp[1], mod.fitted$sp[2], mod.fitted$sp[3], mod.fitted$sp[4]))
        temp.dev.exp<- (deviance(no.temp) - deviance(mod.fitted))/deviance(null.mod)
        return(temp.dev.exp)
      }
    }
    
    dat.train <- dat.train %>%
      mutate(., "mod.fitted.p" = purrr::map(TRAIN.DATA, possibly(m1.p, NA)),
             "mod.fitted.b" = purrr::map(TRAIN.DATA, possibly(m1.b, NA)),
             "temp.devexp.prop.p" = purrr::pmap(list(df = TRAIN.DATA, mod.fitted = mod.fitted.p, response = "PRESENCE"), possibly(temp.dev.exp, NA)),
             "temp.devexp.prop.b" = purrr::pmap(list(df = TRAIN.DATA, mod.fitted = mod.fitted.b, response = "BIOMASS.MOD"), possibly(temp.dev.exp, NA)))
    
    # Merge back the testing list of dataframes
    dat.fitted<- dat.test %>% 
      dplyr::select(., COMNAME, SEASON, TEST.DATA) %>%
      left_join(dplyr::select(dat.train, -TRAIN.TEST), ., by = c("COMNAME", "SEASON"))
    
    # Predictions -------------------------------------------------------------
    predict_func<- function(mod.fitted, response, stat, test.data) {
      temp<- dplyr::select(test.data, one_of(c(pred.vars)))
      test.data<- data.frame(na.omit(temp))
      
      if(response == "PRESENCE") {
        out<- switch(stat,
                     mean = round(as.numeric(predict.gam(mod.fitted, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3),
                     se = round(as.numeric(predict.gam(mod.fitted, newdata = test.data, type = "response", se.fit = TRUE)$se.fit), 3))
        return(out)
      }
      
      if(response == "BIOMASS.MOD") {
        out<- switch(stat,
                     mean = round(exp(as.numeric(predict.gam(mod.fitted, newdata = test.data, type = "response", se.fit = TRUE)$fit)), 3),
                     se = round(exp(as.numeric(predict.gam(mod.fitted, newdata = test.data, type = "response", se.fit = TRUE)$se.fit)), 3))
        return(out)
      }
    }
    
    combined_func<- function(presence, biomass) {
      round(presence*biomass, 3)
    }
    
    # Predict Presence and Biomass
    dat.fitted<- dat.fitted %>%
      mutate("Predicted.mu.p" = purrr::pmap(list(mod.fitted = mod.fitted.p, response = "PRESENCE", stat = "mean", test.data = TEST.DATA), possibly(predict_func, NA)),
             "Predicted.se.p" = purrr::pmap(list(mod.fitted = mod.fitted.p, response = "PRESENCE", stat = "se", test.data = TEST.DATA), possibly(predict_func, NA)),
             "Predicted.mu.b" = purrr::pmap(list(mod.fitted = mod.fitted.b, response = "BIOMASS.MOD", stat = "mean", test.data = TEST.DATA), possibly(predict_func, NA)),
             "Predicted.se.b" = purrr::pmap(list(mod.fitted = mod.fitted.b, response = "BIOMASS.MOD", stat = "se", test.data = TEST.DATA), possibly(predict_func, NA)),
             "Predicted.c" = purrr::map2(Predicted.mu.p, Predicted.mu.b, possibly(combined_func, NA)))
    
    # Calibration and Validation -------------------------------------------------------------
    auc_func<- function(test.data, predicted) {
      library(ROCR)
      temp<- dplyr::select(test.data, one_of(c(pred.vars, "PRESENCE")))
      test.data<- data.frame(na.omit(temp))
      col.ind<- which(colnames(test.data) == "PRESENCE")
      dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
      return(performance(dat, measure = "auc")@y.values[[1]])
    }
    
    dev_exp_func<- function(mod.fitted) {
      return(round(summary(mod.fitted)$dev.expl, 4)*100)
    }
    
    calib_func<- function(test.data, response, predicted) {
      if(response == "PRESENCE"){
        library(ROCR)
        temp<- dplyr::select(test.data, one_of(c(pred.vars, "PRESENCE")))
        test.data<- data.frame(na.omit(temp))
        col.ind<- which(colnames(test.data) == response)
        dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
        return(performance(dat, measure = "rmse")@y.values[[1]])
      }
      if(response == "BIOMASS.MOD") {
        library(hydroGOF)
        test.data<- dplyr::select(test.data, one_of(c(pred.vars, "BIOMASS.MOD")))
        return(rmse(sim = predicted, obs = exp(test.data$BIOMASS.MOD)))
      }
    }
    
    dat.fitted<- dat.fitted %>%
      mutate(., "Dev.Exp.p" = purrr::map(mod.fitted.p, possibly(dev_exp_func, NA)),
             "Dev.Exp.b" = purrr::map(mod.fitted.b, possibly(dev_exp_func, NA)),
             "AUC" = purrr::map2(TEST.DATA, Predicted.mu.p, possibly(auc_func, NA)),
             "RMSE.p" = purrr::pmap(list(test.data = TEST.DATA, response = "PRESENCE", predicted = Predicted.mu.p), possibly(calib_func, NA)),
             "RMSE.b" = purrr::pmap(list(test.data = TEST.DATA, response = "BIOMASS.MOD", predicted = Predicted.c), possibly(calib_func, NA)))
    
    mod.perf.table<- dat.fitted %>%
      dplyr::select(., COMNAME, SEASON, Dev.Exp.p, Dev.Exp.b, temp.devexp.prop.p, temp.devexp.prop.b, AUC, RMSE.p, RMSE.b) %>%
      unnest(.) %>%
      data.frame(.)
    
    # Write em out
    write.csv(mod.perf.table, file = "./Results/GAM.DLN.mod.table.csv")
    
    # Remove data column
    train.dat.col<- match("TRAIN.DATA", colnames(dat.fitted))
    dat.fitted<- dat.fitted[,-train.dat.col]
    
    # Save fitted object
    saveRDS(dat.fitted, file = "./Data/dat.fitted.rds")
    
    rm(dat.fitted, dat.train, dat.test)
    gc()
  }
  
  if(model.framework == "boruta.biomass"){
    boruta_func<- function(df, response) {
      if(response == "PRESENCE") {
        keep.vec<- c(pred.vars, "PRESENCE")
        temp<- dplyr::select(df, one_of(keep.vec))
        mod.dat<- as.data.frame(na.omit(temp))
        return(Boruta(y = as.factor(mod.dat[,"PRESENCE"]), x = mod.dat[,pred.vars], doTrace = 0))
      }
      if(response == "BIOMASS.MOD") {
        keep.vec<- c(pred.vars, "BIOMASS.MOD")
        temp<- dplyr::select(df, one_of(keep.vec))
        mod.dat<- as.data.frame(na.omit(temp))
        return(Boruta(y = mod.dat[,"BIOMASS.MOD"], x = mod.dat[,pred.vars], doTrace = 0))
      }
    }
    
    boruta_vars<- function(boruta.out){
      names(boruta.out$finalDecision[boruta.out$finalDecision %in% c("Confirmed", "Tentative")])
    }
    
    rf_boruta_fit<- function(df, response, boruta.out) {
      if(response == "PRESENCE") {
        ntrees<- 1000
        keep.vec<- c(pred.vars, "PRESENCE")
        temp<- dplyr::select(df, one_of(keep.vec))
        mod.dat<- as.data.frame(na.omit(temp))
        vars.imp <- names(boruta.out$finalDecision[boruta.out$finalDecision %in% c("Confirmed", "Tentative")])
        return(randomForest(y = as.factor(mod.dat[,"PRESENCE"]), x = mod.dat[,vars.imp], ntree = ntrees, importance = TRUE, proximity = TRUE))
      }
      
      if(response == "BIOMASS.MOD"){
        ntrees<- 1000
        keep.vec<- c(pred.vars, "BIOMASS.MOD")
        temp<- dplyr::select(df, one_of(keep.vec))
        mod.dat<- as.data.frame(na.omit(temp))
        vars.imp <- names(boruta.out$finalDecision[boruta.out$finalDecision %in% c("Confirmed", "Tentative")])
        return(randomForest(y = mod.dat[,"BIOMASS.MOD"], x = mod.dat[,vars.imp], ntree = ntrees, importance = TRUE, proximity = TRUE))
      }
    }
    
    dat.train <- dat.train %>%
      mutate(., "boruta.out.p" = purrr::map2(TRAIN.DATA, "PRESENCE", possibly(boruta_func, NA)),
             "boruta.out.b" = purrr::map2(TRAIN.DATA, "BIOMASS.MOD", possibly(boruta_func, NA)),
             "vars.imp.p" = purrr::map(boruta.out.p, possibly(boruta_vars, NA)),
             "vars.imp.b" = purrr::map(boruta.out.b, possibly(boruta_vars, NA)),
             "mod.fitted.p" = purrr::pmap(list(df = TRAIN.DATA, response = "PRESENCE", boruta.out = boruta.out.p), possibly(rf_boruta_fit, NA)),
             "mod.fitted.b" = purrr::pmap(list(df = TRAIN.DATA, response = "BIOMASS.MOD", boruta.out = boruta.out.b), possibly(rf_boruta_fit, NA)))
    
    # Merge back the testing list of dataframes
    dat.fitted<- dat.test %>% 
      dplyr::select(., COMNAME, TEST.DATA) %>%
      left_join(dplyr::select(dat.train, -TRAIN.TEST), ., by = c("COMNAME"))
    
    # Predictions -------------------------------------------------------------
    predict_func<- function(mod.fitted, response, test.data) {
      temp<- dplyr::select(test.data, one_of(c(pred.vars)))
      test.data<- data.frame(na.omit(temp))
      
      if(response == "PRESENCE") {
        return(round(predict(mod.fitted, newdata = test.data, type = "prob")[,2], 3))
      }
      
      if(response == "BIOMASS.MOD") {
        return(round(exp(as.numeric(predict(mod.fitted, newdata = test.data))), 3))
      }
    }
    
    combined_func<- function(presence, biomass) {
      round(presence*biomass, 3)
    }
    
    # Predict Presence and Biomass
    dat.fitted<- dat.fitted %>%
      mutate("Predicted.p" = purrr::pmap(list(mod.fitted = mod.fitted.p, response = "PRESENCE", test.data = TEST.DATA), possibly(predict_func, NA)),
             "Predicted.b" = purrr::pmap(list(mod.fitted = mod.fitted.b, response = "BIOMASS.MOD", test.data = TEST.DATA), possibly(predict_func, NA)),
             "Predicted.c" = purrr::map2(Predicted.p, Predicted.b, possibly(combined_func, NA)))
    
    # Calibration and Validation -------------------------------------------------------------
    auc_func<- function(test.data, predicted) {
      library(ROCR)
      temp<- dplyr::select(test.data, one_of(c(pred.vars, "PRESENCE")))
      test.data<- data.frame(na.omit(temp))
      col.ind<- which(colnames(test.data) == "PRESENCE")
      dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
      return(performance(dat, measure = "auc")@y.values[[1]])
    }
    
    dev_exp_func<- function(mod.fitted, response) {
      if(response == "PRESENCE") {
        return(round(as.numeric(mod.fitted$err.rate[nrow(mod.fitted$err.rate),1])*100, 3))
      }
      if(response == "BIOMASS.MOD"){
        return(round(mod.fitted$rsq[length(mod.fitted$rsq)]*100, 3))
      }
    }
    
    calib_func<- function(test.data, response, predicted) {
      if(response == "PRESENCE"){
        library(ROCR)
        temp<- dplyr::select(test.data, one_of(c(pred.vars, "PRESENCE")))
        test.data<- data.frame(na.omit(temp))
        col.ind<- which(colnames(test.data) == response)
        dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
        return(performance(dat, measure = "rmse")@y.values[[1]])
      }
      if(response == "BIOMASS.MOD") {
        library(hydroGOF)
        test.data<- dplyr::select(test.data, one_of(c(pred.vars, "BIOMASS.MOD")))
        return(rmse(sim = predicted, obs = exp(test.data$BIOMASS.MOD)))
      }
    }
    
    dat.fitted<- dat.fitted %>%
      mutate(., "Dev.Exp.p" = purrr::map2(mod.fitted.p, "PRESENCE", possibly(dev_exp_func, NA)),
             "Dev.Exp.b" = purrr::map2(mod.fitted.b, "BIOMASS.MOD", possibly(dev_exp_func, NA)),
             "AUC" = purrr::map2(TEST.DATA, Predicted.p, possibly(auc_func, NA)),
             "RMSE.p" = purrr::pmap(list(test.data = TEST.DATA, response = "PRESENCE", predicted = Predicted.p), possibly(calib_func, NA)),
             "RMSE.b" = purrr::pmap(list(test.data = TEST.DATA, response = "BIOMASS.MOD", predicted = Predicted.c), possibly(calib_func, NA)))
    
    # Save some results -------------------------------------------------------------
    # Model Validation and Calibration table
    mod.perf.table<- dat.fitted %>%
      dplyr::select(., COMNAME, Dev.Exp.p, Dev.Exp.b, AUC, RMSE.p, RMSE.b) %>%
      unnest(.) %>%
      data.frame(.)
    
    mod.selection.p<- data.frame(matrix(nrow = length(pred.vars), ncol = nrow(dat.fitted)+1))
    names(mod.selection.p)[1]<- "Variable"
    mod.selection.p[,1]<- pred.vars
    
    sv<- 2
    
    for(i in 1:nrow(dat.fitted)) {
      dat.use<- dat.fitted[i,]
      boruta.use<- dat.use$boruta.out.p[[1]]
      
      if(any(!is.na(dat.use$mod.fitted.p[[1]]))) {
        vars.imp<- data.frame(boruta.use["ImpHistory"])
        vars.imp.means<- colMeans(vars.imp, na.rm = T)
        names(vars.imp.means)<- gsub("ImpHistory.", "", names(vars.imp.means))
        col.sv<- as.numeric(na.omit(match(names(vars.imp.means), mod.selection.p$Variable)))
        mod.selection.p[col.sv, sv]<- as.numeric(vars.imp.means[col.sv])
        names(mod.selection.p)[sv]<- paste(simpleCap(as.character(dat.use$COMNAME)), scenario, sep = ".")
        sv<- sv+1
      } else {
        vars.imp<- NA
        mod.selection.p[,sv]<- NA
        names(mod.selection.p)[sv]<- paste(simpleCap(as.character(dat.use$COMNAME)), scenario, sep = ".")
        sv<- sv+1
      }
    }
    
    mod.selection.p$Model<- rep("Presence", nrow(mod.selection.p))
    
    mod.selection.b<- data.frame(matrix(nrow = length(pred.vars), ncol = nrow(dat.fitted)+1))
    names(mod.selection.b)[1]<- "Variable"
    mod.selection.b[,1]<- pred.vars
    
    sv<- 2
    
    for(i in 1:nrow(dat.fitted)) {
      dat.use<- dat.fitted[i,]
      boruta.use<- dat.use$boruta.out.b[[1]]
      
      if(any(!is.na(dat.use$mod.fitted.b[[1]]))) {
        vars.imp<- data.frame(boruta.use["ImpHistory"])
        vars.imp.means<- colMeans(vars.imp, na.rm = T)
        names(vars.imp.means)<- gsub("ImpHistory.", "", names(vars.imp.means))
        col.sv<- as.numeric(na.omit(match(names(vars.imp.means), mod.selection.b$Variable)))
        mod.selection.b[col.sv, sv]<- as.numeric(vars.imp.means[col.sv])
        names(mod.selection.b)[sv]<- paste(simpleCap(as.character(dat.use$COMNAME)), scenario, sep = ".")
        sv<- sv+1
      } else {
        vars.imp<- NA
        mod.selection.b[,sv]<- NA
        names(mod.selection.b)[sv]<- paste(simpleCap(as.character(dat.use$COMNAME)), scenario, sep = ".")
        sv<- sv+1
      }
    }
    
    mod.selection.b$Model<- rep("Biomass", nrow(mod.selection.b))
    mod.selection<- bind_rows(mod.selection.p, mod.selection.b)
    
    # Write em out
    write.csv(mod.perf.table, file = paste(out.dir, scenario, ".mod.table.csv", sep = ""))
    write.csv(mod.selection, file = paste(out.dir, scenario, ".mod.selection.csv", sep = ""))
    
    # Remove data column
    train.dat.col<- match("TRAIN.DATA", colnames(dat.fitted))
    dat.fitted<- dat.fitted[,-train.dat.col]
    
    # Save fitted object
    save(dat.fitted, file = paste(out.dir, scenario, ".dat.fitted.RData", sep = ""))
    
    rm(dat.fitted, dat.train, dat.test)
    gc()
  }
}