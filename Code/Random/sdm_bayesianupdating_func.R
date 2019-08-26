#### SDM Bayes Function: Written by A.Pershing in matlab, converted into R
# Simulate some data using the gamSim function, which we will then fit a model 

SDMbayes<- function(x, bmn, bsd, fmn, fsd, nevaD.neg, nevaD.neut, nevaD.pos, nevaV.low, nevaV.mod, nevaV.high, nevaV.vhigh){
  
  ## Function details
  # x = Vector of prob of occurrence to generate probability desnity function
  # bmn = mean in base case
  # bsd = std. dev of the base case
  # fmn = mean in the forecast
  # fsd = std. dev of the forecast
  # nevaD.neg, nevaD.neut, nevaD.pos = NEVA directional effect voting (-, 0, +)
  # nevaV.low, nevaV.mod, nevaV.high, nevaV.vhigh = NEVA vulnerability voting (1, 2, 3, 4)
  
  if(FALSE){
    x<- seq(from = -5, to = 5, length.out = 500)
    bmn<- 0 #-5
    bsd<- 1
    fmn<- 1 #3
    fsd<- 1
    nevaD.neg = 0
    nevaD.neut = 3
    nevaD.pos = 9
    nevaV.low = 0
    nevaV.mod = 0
    nevaV.high = 3
    nevaV.vhigh = 9
    #nevaD = c(0, 5, 5)
    #nevaV = c(24, 24, 21, 32)
  }
  
  # Preliminaries
  nevaD<- c(nevaD.neg, nevaD.neut, nevaD.pos) # Collect directional effect votes
  nevaV<- c(nevaV.low, nevaV.mod, nevaV.high, nevaV.vhigh) # Collect vulnerability votes
  N<- 100 # Number of "models" to generate -- hacking process, assuming that our observed future is one possible future
  Np<- rep(0, length(x)) # Storage vector -- what is this storing? 
  pT<- rep(0, N) # Storage vector -- what is this storing?
  Dbrks<- c(-1, 1)*bsd+bmn # Directional effect cuts (+/- 1SD), sets bins to integrate future pdf over
  Vbrks<- c(-3, -2, -1, 1, 2, 3)*bsd/2+bmn # Vulnerability cuts, sets bins to integrate future pdf over
  md<- length(Vbrks)/2+1 # Numeric value (2)
  
  # Enter loop to generate "models"
  for(i in 1:N){
    # Create a new model
    # This is the hack. Assuming that the model means vary with the SD of the projection model. Also assuming that the perturbations to the mean increase the variance. Not sure if these make sense for this application, but a way at coming up with the elusive P(model). Note, I wonder if we can use the lpmatrix from GAM model to get this?
    # Generate new mean and new SD for this potential future "model"
    r<- rnorm(1)*fsd # Random pertubation
    fm2<- fmn + r # Perturbation to mean
    fs2<- fsd + abs(r) # Perturbation to SD
    M<- dnorm(x, fm2, fs2) # Our potential future "model". This is our prior -- SDM quantitative projection pdf. 
    
    ## Directional effect piece
    dtmp<- pnorm(Dbrks, fm2, fs2) # Cumulative distribution function, p(x) <= -1 SD and p(x) <= 1SD from N(fm2, fs2) dist
    pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
    
    pD<- dmultinom(nevaD, prob = pd) # Probability of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    
    ## Vulnerability piece
    vtmp<- pnorm(Vbrks, fm2, fs2) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    vd<- rep(NA, md)
    vd[1]<- vtmp[md]
    
    for(k in 1:(md-1)){
      vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      #print(paste("This is ", k, " and vd[k] ", vd[k], "and vd[k+1] ", vd[k+1], sep = ""))
    }
    
    pV<- dmultinom(nevaV, prob = vd) # Probability of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    
    # Combined probability of getting the directional effect votes and the vulnerability votes
    pT[i]<- pD*pV
    
    # Need a bit more explantion here...Np starts at 0 -- we add pT[i] = probability of observing directional effect votes and vulnerability votes, given our future model M, then multiply by M = our potential "future" model (our prior?)
    Np<- Np+pT[i]*M 
  }
  
  totP<- sum(pT)
  
  if(totP<1e-12){
    #print("Total P is pretty much zero")
  }
  
  # Posterior PDF
  Np.out<- Np/totP
  
  # Added step
  #Np.out2<- cNp.out + (1-c)normpdf(x, fmn, fsd)
  
  if(all(is.nan(Np.out))){
    return(NA)
  } else {
    return(Np.out)
  }
}
