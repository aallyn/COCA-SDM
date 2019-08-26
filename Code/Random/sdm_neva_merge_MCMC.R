###### Merging SDM and NEVA: Approach 3.14159
################
## Preliminaries
library(mgcv)
library(tidyverse)
library(mvtnorm)
library(Smisc)
library(rgeos)
library(akima)
library(sf)
library(viridis)
library(cowplot)
library(corrr)


# Step by Step — old ------------------------------------------------------
## Step 1: Data prep. We have a dataset, that includes the point locations of NOAA NEFSC bottom trawl tows, presence/absence of fish, and then a suite of environmental characteristics at these locations including tempetature, depth and along shelf position.
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")
spp.example<- "SUMMER FLOUNDER"

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "BIOMASS.WMEAN.ANOM" = as.numeric(BIOMASS.WMEAN.ANOM),
         "BIOMASS.WMEAN.Scale" = as.numeric(scale(BIOMASS.WMEAN)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST)),
         "BOTTOMTEMP.Scale" = as.numeric(scale(BOTTEMP))) %>%
  filter(., COMNAME == "ATLANTIC COD")

train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither")) 

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


## Step 2: The species distribution modeling. Use GAM to correlate environmental characteristics with presence/absence data.
gam.mod<- gam(PRESENCE ~ s(SHELF_POS.Scale, fx = FALSE, bs = 'cs') + s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train$TRAIN.DATA[[2]], family = binomial(link = logit), select = TRUE, weights = ABUNDANCE)

## Step 3: Investigating the details of the fitted GAM. 
summary(gam.mod)
plot(gam.mod, pages = 1)

coef.mod<- coef(gam.mod) # Coefficients for every term/basis function in the model
coef.mod
vc.mod<- vcov(gam.mod) # Variance/covariance matrix
vc.mod
sqrt(diag(vc.mod)) # The standard errors

## Step 3: Simulating potential outcomes from the model, rather than just the mean
# For a LM model, the point estimates and the SE describe a multivariate normal distribution (could this get us the likelihood we are after??). We can then imagine pulling samples from that multivariate distribution (each pull would be a different parameter estimate), and then use those coefficients to make our predictions. For the GAM, we are after something similar, but things are a bit more complicated with the basis functions. So, instead of multiplying each simulation, we generate the linear predictor matrix for the observations and multiply that by the model coefficients to get simulations. 
lp.mod<- predict(gam.mod, type = "lpmatrix") # Note, here we could supply new data
str(lp.mod) # Dimensions: nrow = number of points, ncol = a column for every basis function, plus the constant term

# Simulating from a multivariate normal to get different "models" (different basis functions/smooths)
set.seed(35)
mod.sim <- rmvnorm(10000, mean = coef.mod, sigma = vc.mod) # Dimensions: nrow = number of draws, ncol = a column for every basis function, plus the constant term  

mod.fits<- lp.mod %*% t(mod.sim) # Fitted values - this is what we would synthesize across the shelf and then confront with the NEVA data

## Step 3.5: We can definitely generate a number (any number) of potential "models". Is there a way that we can calculate the likelihood of each of these models?
# Example
mod.sim1<- rmvnorm(1, mean = coef.mod, sigma = vc.mod) # Random sample....
pmod.sim1<- dmvnorm(x = mod.sim1, mean = coef.mod, sigma = vc.mod, log = TRUE) # I think this is doing what we want...

pmod.sim<- dmvnorm(x = mod.sim, mean = coef.mod, sigma = vc.mod, log = TRUE)

plot.dat<- data.frame("Model" = seq(from = 1, to = length(pmod.sim), by = 1), "Relative.Log.Likelihood" = pmod.sim)
plot.dat<- dplyr::arrange(plot.dat, Relative.Log.Likelihood)
plot.dat$Plot.Ind<- seq(from = 1, to = nrow(plot.dat))

plot.out<- ggplot(plot.dat, aes(x = Plot.Ind, y = Relative.Log.Likelihood)) +
  geom_line(color = "blue") +
  theme_bw()
plot.out

# Want the model that maximizes the log likelihood?
mod.maxll.ind<- plot.dat$Model[which.max(plot.dat$Relative.Log.Likelihood)]
mod.maxll<- mod.sim[mod.maxll.ind,]
mod.maxll<- matrix(mod.maxll, nrow = 1, ncol = length(coef.mod), byrow = T, dimnames = list(NULL, names(coef(gam.mod))))

# Predictions from this model
mod.fit.best<- lp.mod %*% t(mod.maxll) 



# M-H MCMC with our data --------------------------------------------------
## M-H Algorithm for our example
# NegativeLL function -- the likelihood of NEVA votes given the base and future, which comes out of a potential gam model -- combined
loglike_func<- function(gam.mod, season, cand.params, base.preds, fut.preds, nevaD, nevaV){
  
  if(FALSE){
    gam.mod = mod
    season = "Fall"
    cand.params = rmvnorm(1, mean = coef.mod, sigma = vc.mod) 
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = c(0, 3, 9)
    nevaV = c(0, 0, 3, 9)
  }
  
  # NEVA conversion components
  # New baseline and future for specific candidate parameters
  # Get linear predictor matrix for each of the observations in base.preds and fut.preds
  base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
  fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  lpmat.base<- predict.gam(gam.mod, newdata = base.preds.use, type = "lpmatrix")
  lpmat.fut<- predict.gam(gam.mod, newdata = fut.preds.use, type = "lpmatrix")
  
  # Predictions, combining lpmats and candidate parameters
  # Multiply candidates by linear predictor matrix
  base.prob<- lpmat.base %*% coef(gam.mod)
  fut.prob<- lpmat.fut %*% t(cand.params)
  
  # For proposed GAM, get new baseline and future statistics
  bmn<- mean(base.prob, na.rm = T)
  bsd<- sd(base.prob, na.rm = T)
  fmn<- mean(fut.prob, na.rm = T)
  fsd<- sd(fut.prob, na.rm = T)
  
  # Bin baseline distribution for directional effects and then bin vulnerability
  Dbrks<- c(-1, 1)*bsd+bmn # Directional effect cuts (+/- 1SD), sets bins to integrate future pdf over
  Vbrks<- c(-3, -2, -1, 1, 2, 3)*bsd/2+bmn # Vulnerability cuts, sets bins to integrate future pdf over
  md<- length(Vbrks)/2+1 # Numeric value (2)
  
  ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
  if(!is.null(nevaD)){
    dtmp<- pnorm(Dbrks, fmn, fsd) # Cumulative distribution function, p(x) <= -1 SD and p(x) <= 1SD from N(fm2, fs2) dist
    pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
    
    like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
  }

  ## Vulnerability piece
  if(!is.null(nevaV)){
    vtmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    vd<- rep(NA, md)
    vd[1]<- vtmp[md]
    
    for(k in 1:(md-1)){
      vd[k+1]<- vtmp[md+k]+vtmp[md-k]
    }
    
    like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
  }
  
  if(is.null(nevaD)){
    likelihood<- like.vuln
  } else if(is.null(nevaV)){
    likelihood<- like.dir
  } else {
    likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
  }
  return(likelihood) # Loglikelihood
}

# Probability of our model 
logprior_func<- function(gam.mod, cand.params){
  # Likelihood of the GAM baseline future (fitted curves)
  # Gather coefficiences and varcov matrix, needed to define the normal dist of each candidate parameter value
  coefs.mod<- coef(gam.mod) # Coefficients for every term/basis function in the model
  vcov.mod<- vcov(gam.mod) 
  
  # Likelihood of these candidate parameter values, given coef (mean) and sigma (se)
  cand.params.logprior<- dmvnorm(cand.params, mean = coefs.mod, sigma = vcov.mod, log = TRUE)
  return(cand.params.logprior)
}

# Posterior (loglike + logprior)
logposterior_func<- function(gam.mod, season, cand.params, base.preds, fut.preds, nevaD, nevaV){
  # Likelihood of NEVA dir and vuln votes, given GAM baseline future (fitted curves)
  if(class(cand.params) == "numeric"){
    cand.params<- matrix(cand.params, nrow = 1, ncol = length(cand.params), byrow = T, dimnames = list(NULL, names(coef(gam.mod))))
  }
  
  loglike.out<- loglike_func(gam.mod, season, cand.params, base.preds, fut.preds, nevaD, nevaV)
  
  # Probability of prior
  logprior.out<- logprior_func(gam.mod, cand.params)
  
  # Add em all together
  return(data.frame("Likelihood" = loglike.out, "Prior" = logprior.out, "Posterior" = loglike.out + logprior.out))
}

# Proposal function (how do we get our new candidate parameters)
proposal_func<- function(cand.params, gam.mod, tune){
  # Each one independently -- not working well
  # prop.out<- t(as.matrix(rnorm(length(cand.params), mean = cand.params, sd = rep(tune, length(cand.params)))))
  
  # Multivariate normal -- works, but odd
  prop.out<- rmvnorm(1, mean = cand.params, sigma = vcov(gam.mod)*tune)
  rownames(prop.out)<- NULL
  colnames(prop.out)<- names(coef(gam.mod))
  return(prop.out)
  
  # Something else...https://stats.stackexchange.com/questions/91682/how-to-draw-estimates-based-on-variance-covariance-matrix
  prop.out<- t(as.matrix(rnorm(length(cand.params), mean = cand.params, tune) + chol(vcov(gam.mod))%*%rnorm(length(cand.params), tune)))
  return(prop.out)
}

# MCMC MH Algorithm
mh_mcmc<- function(gam.mod, season, cand.params.start, base.preds, fut.preds, nevaD, nevaV, iterations, tune){
  
  out.likes<- array(dim = c(iterations, 3))
  out.cands<- array(dim = c(iterations+1, length(cand.params.start)))

  out.cands[1,]<- cand.params.start
  
  for(i in 1:iterations){
    
    # Print progress
    if(i %% 1000 == 0){
      cat("iter", i, "\n")
    }
    
    # Update candidate parameter values
    proposal<- proposal_func(out.cands[i,], gam.mod, tune)
    
    # Calculate likelihood, prior and posterior given proposal and previous candidate
    post.curr<- logposterior_func(gam.mod, season, cand.params = proposal, base.preds, fut.preds, nevaD, nevaV)
    post.prev<- logposterior_func(gam.mod, season, cand.params = out.cands[i,], base.preds, fut.preds, nevaD, nevaV)
    
    # Store likelihood and prior
    out.likes[i,]<- cbind(post.curr$Likelihood, post.curr$Prior, post.curr$Posterior)
    
    # Compare current and previous
    probab<- exp(post.curr$Posterior - post.prev$Posterior)
    
    # Keep new proposal or previous one, with some randomness added to "jump"
    if(runif(1) < probab){
      out.cands[i+1,]<- proposal
    } else{
      out.cands[i+1,]<- out.cands[i,]
    }
  }
  
  mh.mcmc.result<- list(out.likes, out.cands)
  return(mh.mcmc.result)
}

# Running the mh_mcmc algorithm. Usually want to do all of the process with two chains starting from two different points to ensure that the algorithm converges.
# Need base.preds, fut.preds, nevaD and nevaV
spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds03232018.rds"
spring.preds<- readRDS(spring.preds)

fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds03232018.rds"
fall.preds<- readRDS(fall.preds)

base.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Baseline)),
         "SEASON" = rep("SPRING", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Baseline)),
         "SEASON" = rep("FALL", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

fut.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Spring.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Spring.2055)),
         "SEASON" = rep("SPRING", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Fall.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Fall.2055)),
         "SEASON" = rep("FALL", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

#### A few examples
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")
spp.example<- c("ATLANTIC COD", "SUMMER FLOUNDER", "LONGFIN SQUID", "AMERICAN LOBSTER")

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "BIOMASS.WMEAN.ANOM" = as.numeric(BIOMASS.WMEAN.ANOM),
         "BIOMASS.WMEAN.Scale" = as.numeric(scale(BIOMASS.WMEAN)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST)),
         "BOTTOMTEMP.Scale" = as.numeric(scale(BOTTEMP))) %>%
  filter(., COMNAME %in% spp.example)
dat$PRESENCE.ABUNDANCE<- ifelse(dat$ABUNDANCE > 0, 1, 0)
dat$WT.ABUNDANCE<- ifelse(dat$PRESENCE.ABUNDANCE == 0, 1, dat$ABUNDANCE)

train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither")) 

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

## Loop over a few different scenarios...
dir.scenarios<- data.frame("Dir.Scenario" = c("Negative", "Neutral", "Positive"), "Negative" = c(100, 10, 0), "Neutral" = c(20, 100, 20), "Positive" = c(0, 10, 100))
nevaV<- NULL

vuln.scenarios<- data.frame("Vuln.Scenario" = c("Low", "Mod", "High", "VeryHigh"), "Low" = c(10, 1, 0, 0), "Mod" = c(2, 10, 1, 0), "High" = c(0, 1, 10, 2), "VeryHigh" = c(0, 0, 1, 10))
nevaD<- NULL

dir.scenarios.c<- data.frame("Dir.Scenario" = c("Negative", "Positive"), "Negative" = c(10, 0), "Neutral" = c(2, 2), "Positive" = c(0, 10))
vuln.scenarios.c<- data.frame("Vuln.Scenario" = c("Low", "VeryHigh"), "Low" = c(10, 0), "Mod" = c(2, 0), "High" = c(0, 2), "VeryHigh" = c(0, 10))
scenarios.c<- expand.grid(dir.scenarios.c$Dir.Scenario, vuln.scenarios.c$Vuln.Scenario)
names(scenarios.c)<- c("Dir.Scenario", "Vuln.Scenario")

scenarios.c<- scenarios.c %>%
  left_join(., dir.scenarios.c, by = "Dir.Scenario")

scenarios.c<- scenarios.c %>%
  left_join(., vuln.scenarios.c, by = "Vuln.Scenario")

scenarios.c$Scenario<- paste(scenarios.c$Dir.Scenario, scenarios.c$Vuln.Scenario, sep = ".")

n.chains<- 2
burn.in<- 3000
iterations<- 7000
thin<- 3
set.seed(13)

scenarios.run<- dir.scenarios

out.dir.stem<- "~/Desktop/MCMC_Results5/"

## Loop
for(i in 1:nrow(dat.train)){
  dat.use<- dat.train[i,]
  season.use<- as.character(dat.use$SEASON[[1]])
  gam.mod<- gam(PRESENCE.ABUNDANCE ~ s(SHELF_POS.Scale, fx = FALSE, k = 5, bs = 'cs') + s(DEPTH.Scale, k = 5, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, k = 5, fx = FALSE, bs = 'cs'), weights = WT.ABUNDANCE, drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)
  coef.mod<- coef(gam.mod)
  vc.mod<- vcov(gam.mod)
  
  for(j in 1:nrow(scenarios.run)){
    scenario.use<- scenarios.run[j,]
    cand.params.start<- list(rmvnorm(1, mean = rnorm(n = length(coef.mod), mean = coef.mod, sd = 3), sigma = vc.mod),
                             rmvnorm(1, mean = rnorm(n = length(coef.mod), mean = coef.mod, sd = 3), sigma = vc.mod))
    mcmc.results<- vector("list", n.chains)
    
    for(k in 1:n.chains){
      cand.params.use<- cand.params.start[[k]]
      rownames(cand.params.use)<- NULL
      colnames(cand.params.use)<- names(coef.mod)
      
      mcmc.results[[k]]<- mh_mcmc(gam.mod = gam.mod, season = season.use, cand.params.start = cand.params.use, base.preds = base.preds, fut.preds = fut.preds, nevaD = c(scenario.use$Negative, scenario.use$Neutral, scenario.use$Positive), nevaV = NULL, iterations = iterations, tune = 0.6)
    }
    
    # Acceptance rate
    accept.rate<- 1-mean(duplicated(mcmc.results[[1]][2][[1]][-(1:burn.in),]))
    print(paste(tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), " accept rate = ", accept.rate, sep = ""))
    
    # Processing and results plots
    likes<- list(mcmc.results[[1]][[1]], mcmc.results[[2]][[1]])
    likes.thin<- lapply(likes, FUN = function(x) data.frame(x))
    likes.thin<- lapply(likes.thin, FUN = function(x) x[seq(1, nrow(x), thin),])
    likes.nsamps<- nrow(likes.thin[[1]])
    
    posts<- list(mcmc.results[[1]][[2]][burn.in:iterations,], mcmc.results[[2]][[2]][burn.in:iterations,])
    posts.thin<- lapply(posts, FUN = function(x) data.frame(x))
    posts.thin<- lapply(posts.thin, FUN = function(x) x[seq(1, nrow(x), thin),])
    posts.nsamps<- nrow(posts.thin[[1]])
    
    # Likes dataframe and plots of mixing
    likes.df<- data.frame(do.call(rbind, likes.thin))
    names(likes.df)<- c("Likelihood (NEVA|SDM)", "Prior P(SDM)", "Posterior")
    likes.df$Iteration<- rep(seq(from = 1, to = likes.nsamps, by = 1), 2)
    likes.df$Chain<- rep(c("Chain1", "Chain2"), each = likes.nsamps)
    likes.df<- likes.df %>%
      gather(., Sample, Value, -Iteration, -Chain)
    likes.df$Sample<- factor(likes.df$Sample, levels = c("Likelihood (NEVA|SDM)", "Prior P(SDM)", "Posterior"))
    
    out.plot<- ggplot(likes.df) +
      geom_line(aes(x = Iteration, y = Value, color = Chain), alpha = 0.25) +
      #scale_color_manual(values = c('#4daf4a', '#377eb8', '#e41a1c')) + 
      ylab("Log Likelihood") + 
      facet_wrap(~Sample, scales = "free") + 
      theme_bw()
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_LikePriorPost.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
    
    # Plots
    plot.dat<- data.frame(do.call(rbind, posts.thin))
    colnames(plot.dat)<- names(coef(gam.mod))
    plot.dat$Iteration<- rep(seq(from = 1, to = posts.nsamps, by = 1), 2)
    plot.dat$Chain<- rep(c("Chain1", "Chain2"), each = posts.nsamps)
    plot.dat<- plot.dat %>%
      gather(., Parameter, Value, -Iteration, -Chain)
    plot.dat$Parameter<- factor(plot.dat$Parameter, levels = names(coef(gam.mod)))
    plot.dat$Chain<- factor(plot.dat$Chain, levels = c("Chain1", "Chain2"))
    
    out.plot<- ggplot(plot.dat) + 
      geom_line(aes(x = Value, color = Chain), stat = "density") +
      facet_wrap(~Parameter, nrow = 4, scales = "free") + 
      theme_bw()
    
    # Add original means
    orig<- data.frame("Parameter" = names(coef.mod), "Value" = coef.mod)
    
    out.plot2<- out.plot +
      geom_vline(data = orig, aes(xintercept = Value), color = "black") +
      facet_wrap(~Parameter, nrow = 4, scales = "free") 
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_PostDist.jpg", sep = ""), out.plot2, width = 11, height = 8, dpi = 125, units = "in")
    
    # Iteration series...
    out.plot<- ggplot(plot.dat) + 
      geom_line(aes(x = Iteration, y = Value, color = Chain, group = Chain), alpha = 0.4) +
      scale_x_continuous(breaks = seq(from = 1, to = iterations, by = 1000)) +
      facet_wrap(~Parameter, nrow = 4, scales = "free") + 
      theme_bw()
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_PostMix.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")

    # Maps
    source("~/GitHub/COCA/Code/plot_func.R")
    maps<- plot_func(gam.mod = gam.mod, season = season.use, base.preds = base.preds, fut.preds = fut.preds, out = posts.thin)
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_Maps.jpg", sep = ""), maps, width = 11, height = 8, dpi = 125, units = "in")
    rm(maps)
    
    # Save MCMC results
    saveRDS(mcmc.results, file = paste(out.dir.stem, "mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), ".rds", sep = ""))
  }
}


### Parallel processing?
# Need these results faster, run using multiple cores with a core per species?
mcmc_func_wrapper<- function(i, dat.train, scenarios.run = scenarios.run, n.chains = 2, burn.in = 2500, iterations = 10000, thin = 3) {
  
  dat.use<- dat.train[i,]
  season.use<- as.character(dat.use$SEASON[[1]])
  gam.mod<- gam(PRESENCE.ABUNDANCE ~ s(DEPTH.Scale, k = 5, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, k = 5, fx = FALSE, bs = 'cs'), weights = WT.ABUNDANCE, drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)
  coef.mod<- coef(gam.mod)
  vc.mod<- vcov(gam.mod)
  
  for(j in 1:nrow(scenarios.run)){
    scenario.use<- scenarios.run[j,]
    cand.params.start<- list(rmvnorm(1, mean = rnorm(n = length(coef.mod), mean = coef.mod, sd = 3), sigma = vc.mod),
                             rmvnorm(1, mean = rnorm(n = length(coef.mod), mean = coef.mod, sd = 3), sigma = vc.mod))
    mcmc.results<- vector("list", n.chains)
    
    for(k in 1:n.chains){
      cand.params.use<- cand.params.start[[k]]
      rownames(cand.params.use)<- NULL
      colnames(cand.params.use)<- names(coef.mod)
      
      mcmc.results[[k]]<- mh_mcmc(gam.mod = gam.mod, season = season.use, cand.params.start = cand.params.use, base.preds = base.preds, fut.preds = fut.preds, nevaD = c(scenario.use$Negative, scenario.use$Neutral, scenario.use$Positive), nevaV = NULL, iterations = iterations, tune = 0.6)
    }
    
    # Acceptance rate
    accept.rate<- 1-mean(duplicated(mcmc.results[[1]][2][[1]][-(1:burn.in),]))
    print(paste(tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), " accept rate = ", accept.rate, sep = ""))
    
    # Processing and results plots
    likes<- list(mcmc.results[[1]][[1]], mcmc.results[[2]][[1]])
    likes.thin<- lapply(likes, FUN = function(x) data.frame(x))
    likes.thin<- lapply(likes.thin, FUN = function(x) x[seq(1, nrow(x), thin),])
    likes.nsamps<- nrow(likes.thin[[1]])
    
    posts<- list(mcmc.results[[1]][[2]][burn.in:iterations,], mcmc.results[[2]][[2]][burn.in:iterations,])
    posts.thin<- lapply(posts, FUN = function(x) data.frame(x))
    posts.thin<- lapply(posts.thin, FUN = function(x) x[seq(1, nrow(x), thin),])
    posts.nsamps<- nrow(posts.thin[[1]])
    
    # Likes dataframe and plots of mixing
    likes.df<- data.frame(do.call(rbind, likes.thin))
    names(likes.df)<- c("Likelihood (NEVA|SDM)", "Prior P(SDM)", "Posterior")
    likes.df$Iteration<- rep(seq(from = 1, to = likes.nsamps, by = 1), 2)
    likes.df$Chain<- rep(c("Chain1", "Chain2"), each = likes.nsamps)
    likes.df<- likes.df %>%
      gather(., Sample, Value, -Iteration, -Chain)
    likes.df$Sample<- factor(likes.df$Sample, levels = c("Likelihood (NEVA|SDM)", "Prior P(SDM)", "Posterior"))
    
    out.plot<- ggplot(likes.df) +
      geom_line(aes(x = Iteration, y = Value, color = Chain), alpha = 0.25) +
      #scale_color_manual(values = c('#4daf4a', '#377eb8', '#e41a1c')) + 
      ylab("Log Likelihood") + 
      facet_wrap(~Sample, scales = "free") + 
      theme_bw()
    ggsave(paste("~/Desktop/MCMC_Results3/", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_LikePriorPost.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
    
    # Plots
    plot.dat<- data.frame(do.call(rbind, posts.thin))
    colnames(plot.dat)<- names(coef(gam.mod))
    plot.dat$Iteration<- rep(seq(from = 1, to = posts.nsamps, by = 1), 2)
    plot.dat$Chain<- rep(c("Chain1", "Chain2"), each = posts.nsamps)
    plot.dat<- plot.dat %>%
      gather(., Parameter, Value, -Iteration, -Chain)
    plot.dat$Parameter<- factor(plot.dat$Parameter, levels = names(coef(gam.mod)))
    plot.dat$Chain<- factor(plot.dat$Chain, levels = c("Chain1", "Chain2"))
    
    out.plot<- ggplot(plot.dat) + 
      geom_line(aes(x = Value, color = Chain), stat = "density") +
      facet_wrap(~Parameter, nrow = 4, scales = "free") + 
      theme_bw()
    
    # Add original means
    orig<- data.frame("Parameter" = names(coef.mod), "Value" = coef.mod)
    
    out.plot2<- out.plot +
      geom_vline(data = orig, aes(xintercept = Value), color = "black") +
      facet_wrap(~Parameter, nrow = 4, scales = "free") 
    ggsave(paste("~/Desktop/MCMC_Results3/", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_PostDist.jpg", sep = ""), out.plot2, width = 11, height = 8, dpi = 125, units = "in")
    
    # Iteration series...
    out.plot<- ggplot(plot.dat) + 
      geom_line(aes(x = Iteration, y = Value, color = Chain, group = Chain), alpha = 0.4) +
      scale_x_continuous(breaks = seq(from = 1, to = iterations, by = 1000)) +
      facet_wrap(~Parameter, nrow = 4, scales = "free") + 
      theme_bw()
    ggsave(paste("~/Desktop/MCMC_Results2/", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_PostMix.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
    
    # Maps
    source("~/GitHub/COCA/Code/plot_func.R")
    maps<- plot_func(gam.mod = gam.mod, season = season.use, base.preds = base.preds, fut.preds = fut.preds, out = posts.thin)
    ggsave(paste("~/Desktop/MCMC_Results3/", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_Maps.jpg", sep = ""), maps, width = 11, height = 8, dpi = 125, units = "in")
    rm(maps)
    
    # Save MCMC results
    saveRDS(mcmc.results, file = paste("~/Desktop/MCMC_Results3/mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), ".rds", sep = ""))
  }
}

## Set up cluster
library(foreach)
library(doParallel)
library(snow)
library(doSNOW)
library(tcltk)

# Set up cores
cores.avail<- detectCores() 

# get the cluster going
cl <- makeSOCKcluster(cores.avail-2)
registerDoSNOW(cl)

# Progress bar stuff
pb<- tkProgressBar(max = nrow(dat.train))
progress<- function(n) setTkProgressBar(pb, n)
opts<- list(progress = progress)

# Run sdw_func_wrapper function in parallel
mcmc_runs<- foreach(i = 1:nrow(dat.train), .packages = c("tidyverse", "mgcv", "mvtnorm", "Smisc", "rgeos", "akima", "sf", "viridis", "cowplot", "corrr", "tcltk"), .combine = rbind.data.frame, .options.snow = opts) %dopar%
{
  mcmc_func_wrapper(i, dat.train, scenarios.run = scenarios.run, n.chains = 2, burn.in = 2500, iterations = 10000, thin = 3)
}


### Fitted GAM smooth function plots...
t<- list.files(path = "~/Desktop/MCMC_Results4/", pattern = ".rds")

# Prediction dataframe
pred.dat.f<- with(fut.preds$Data[[1]],
     data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("FALL", 1500)))
rescaled.dat.f<- with(fut.preds$Data[[1]],
                      data.frame("SST" = seq(min(Fall.2055, na.rm = T), max(Fall.2055, na.rm = T), length = 500),
                                 "Depth" = seq(min(DEPTH, na.rm = T), max(DEPTH, na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(SHELF_POS, na.rm = T), max(SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("FALL", 500)))

pred.dat.s<- with(fut.preds$Data[[2]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("SPRING", 1500)))
rescaled.dat.s<- with(fut.preds$Data[[2]],
                      data.frame("SST" = seq(min(Spring.2055, na.rm = T), max(Spring.2055, na.rm = T), length = 500),
                                 "Depth" = seq(min(DEPTH, na.rm = T), max(DEPTH, na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(SHELF_POS, na.rm = T), max(SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("SPRING", 500)))

pred.dat<- bind_rows(pred.dat.f, pred.dat.s)
rescaled.dat<- bind_rows(rescaled.dat.f, rescaled.dat.s)

## Original predictions...
gam_fit_func<- function(mod.data){
  mod<- gam(PRESENCE ~ s(SHELF_POS.Scale, k = 5, fx = FALSE, bs = 'cs') + s(DEPTH.Scale, fx = FALSE, k = 5, bs = 'cs') + s(SEASONALMU.OISST.Scale, k = 5, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = mod.data, family = binomial(link = logit), select = TRUE)
  return(mod)
}

gam_pred_func<- function(mod, preds){
  preds<- predict(mod, preds, type = "response")
  return(preds)
}

dat.train.f<- dat.train %>%
  dplyr::filter(., SEASON == "FALL") %>%
  mutate(., Mod = map(TRAIN.DATA, gam_fit_func),
         Pred = map2(Mod, list(pred.dat.f), gam_pred_func))
dat.train.s<- dat.train %>%
  dplyr::filter(., SEASON == "SPRING") %>%
  mutate(., Mod = map(TRAIN.DATA, gam_fit_func),
         Pred = map2(Mod, list(pred.dat.s), gam_pred_func))

orig<- bind_rows(dat.train.f, dat.train.s) %>%
  mutate(., "Spp.Season" = paste(tolower(COMNAME), tolower(SEASON), sep = "_"))

Mode<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

spp.seasons<- paste(rep(spp.example, each = 2), "_", c("FALL", "SPRING"), sep = "")

for(h in seq_along(spp.seasons)){
  spp.season.use<- tolower(spp.seasons[h])
  mod.use<- orig[match(spp.season.use, orig$Spp.Season),]
  files.list<- list.files(path = "~/Desktop/MCMC_Results4/", pattern = ".rds")[grep(spp.season.use, list.files(path = "~/Desktop/MCMC_Results4/", pattern = ".rds"))]

  for(i in seq_along(files.list)){
    dat.use<- readRDS(paste("~/Desktop/MCMC_Results4/", files.list[[i]], sep = ""))
    
    posts<- list(dat.use[[1]][2][[1]][burn.in:iterations,], dat.use[[2]][2][[1]][burn.in:iterations,])
    posts.thin<- lapply(posts, FUN = function(x) data.frame(x))
    posts.thin<- lapply(posts, FUN = function(x) x[seq(1, nrow(x), thin),])
    
    best.fit<- apply(posts.thin[[1]], 2, Mode)
    best.fit<- matrix(best.fit, nrow = 1, ncol = length(best.fit), byrow = T, dimnames = list(NULL, names(coef(mod.use$Mod[[1]]))))
    
    # Make predictions with these values
    newdata.use<- pred.dat[pred.dat$SEASON == mod.use$SEASON[[1]],]
    ilink <- family(mod.use$Mod[[1]])$linkinv
    lpmat.fut<- predict.gam(mod.use$Mod[[1]], newdata = newdata.use, type = "lpmatrix")
    pred.mcmc<- ilink(lpmat.fut %*% t(best.fit))
    
    orig.fit<- as.numeric(coef(mod.use$Mod[[1]]))
    orig.fit<- matrix(orig.fit, nrow = 1, ncol = length(orig.fit), byrow = T, dimnames = list(NULL, names(coef(mod.use$Mod[[1]]))))
    pred.orig<-  ilink(lpmat.fut %*% t(orig.fit))
    
    rescaled.dat.use<- rescaled.dat[rescaled.dat$Season == mod.use$SEASON[[1]],]

    # Plotting
    want<- 1:500
    sst.dat<- data.frame("Parameter" = rep("SST", 1000), "Value" = rep(rescaled.dat.use$SST, 2), "Pred" = c(pred.orig[want], pred.mcmc[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
    
    want<- 501:1000
    depth.dat<- data.frame("Parameter" = rep("Depth", 1000), "Value" = rep(rescaled.dat.use$Depth, 2), "Pred" = c(pred.orig[want], pred.mcmc[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
    
    want<- 1001:1500
    shelf.dat<- data.frame("Parameter" = rep("ShelfPos", 1000), "Value" = rep(rescaled.dat.use$Shelf_Pos, 2), "Pred" = c(pred.orig[want], pred.mcmc[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
    
    dat.all<- bind_rows(sst.dat, depth.dat, shelf.dat)
    dat.all$Parameter<- factor(dat.all$Parameter, levels = c("SST", "Depth", "ShelfPos"))
    
    out.plot<- ggplot(dat.all, aes(x = Value, y = Pred, group = Model)) + 
      geom_line(aes(color = Model), alpha = 0.75, show.legend = F) +
      ylim(c(0,1)) +
      scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
      theme_bw() +
      facet_wrap(~Parameter, scales = "free")
    
    ggsave(paste("~/Desktop/MCMC_Results4/", gsub("mcmc_", "", gsub(".rds", "", files.list[[i]])), "_PredCurv.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
  }
}















## Differences in smooths


# GAM model coefficients (original) and then best fit coefficients (from MCMC) multiplied by pdat...
pred.orig<- predict(gam.mod, pdat, type = "response")

# Get maximimum values from each parameter estimate
Mode<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}






































nevaD<- c(0, 3, 9)
nevaV<- c(0, 0, 3, 9)

# MH Algorithm
niter<- 50000
out<- matrix(NA, niter, length(coef.mod), dimnames = list(NULL, names(coef.mod)))

# Initialize parameters
cand.params<- rmvnorm(1, mean = coef.mod, sigma = vc.mod) # Random sample....

# Current values of log(posterior)
logpost.curr<- log.posterior(gam.mod = mod, cand.params = cand.params.init, base.preds = base.preds, fut.preds = fut.preds, nevaD = nevaD, nevaV = nevaV)

# Run MCMC algorithm
for(i in 1:niter){
  # Report progress
  if(i %% 1000 == 0)
    cat("iter", i, "\n")
  
  # Update candidate values
  cand.params<- rmvnorm(1, cand.params, 0.1*diag(length(cand.params))) # 0.3 is a tuning parameter
  
  # Evaluate the log(posterior)
  logpost.cand<- log.posterior(gam.mod = mod, cand.params = cand.params, base.preds = base.preds, fut.preds = fut.preds, nevaD = nevaD, nevaV = nevaV)
  
  # Compute Metropolis acceptance probability, r
  r<- exp(logpost.cand - logpost.curr)
  
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    cand.params<- cand.params
    logpost.curr<- logpost.cand
  }
  
  out[i,]<- cand.params # Save samples from each iteration
}

# Plots
plot.dat<- out %>%
  data.frame(.) 
colnames(plot.dat)<- names(coef(mod))
plot.dat<- plot.dat %>%
  gather(., Parameter, Value)
plot.dat$Parameter<- factor(plot.dat$Parameter, levels = names(coef(mod)))

out.plot<- ggplot(plot.dat, aes(Value, group = Parameter)) + 
  geom_histogram(aes(y = ..density..)) +
  facet_wrap(~Parameter, nrow = 4, scales = "free") + 
  theme_bw()
out.plot

# Add curve
normdens_func<- function(df){
  grid<- seq(min(df$Value), max(df$Value), length = 10000)
  normdens<- dnorm(grid, mean(df$Value), sd(df$Value))
  return(data.frame(Value = grid, Density = normdens))
}

normdens.plot<- plot.dat %>%
  group_by(., Parameter) %>%
  nest() %>%
  mutate(., "NormDens" = purrr::map(data, normdens_func)) %>%
  dplyr::select(., Parameter, NormDens) %>%
  unnest() %>%
  data.frame()

out.plot2<- out.plot + 
  geom_line(aes(y = Density), data = normdens.plot, color = "red") +
  facet_wrap(~Parameter, nrow = 4, scales = "free") 

# Add original means
orig<- data.frame("Parameter" = names(coef.mod), "Value" = coef.mod)

out.plot3<- out.plot2 +
  geom_vline(data = orig, aes(xintercept = Value), color = "blue") +
  facet_wrap(~Parameter, nrow = 4, scales = "free") 

## Maps
source("./Code/plot_func.R")
maps<- plot_func(gam.mod = gam.mod, base.preds = base.preds, fut.preds = fut.preds, out = chain)









