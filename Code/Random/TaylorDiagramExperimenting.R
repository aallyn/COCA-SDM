dat <- selectByDate(mydata, year = 2003)
dat <- data.frame(date = mydata$date, obs = mydata$nox, mod = mydata$nox)

## now make mod worse by adding bias and noise according to the month
## do this for 3 different models
dat <- transform(dat, month = as.numeric(format(date, "%m")))
mod1 <- transform(dat, mod = mod + 10 * month + 10 * month * rnorm(nrow(dat)),
                  model = "model 1")
## lag the results for mod1 to make the correlation coefficient worse
## without affecting the sd
mod1 <- transform(mod1, mod = c(mod[5:length(mod)], mod[(length(mod) - 3) :
                                                          length(mod)]))

## model 2
mod2 <- transform(dat, mod = mod + 7 * month + 7 * month * rnorm(nrow(dat)),
                  model = "model 2")
## model 3
mod3 <- transform(dat, mod = mod + 3 * month + 3 * month * rnorm(nrow(dat)),
                  model = "model 3")

mod.dat <- rbind(mod1, mod2, mod3)

## basic Taylor plot
TaylorDiagram(mod.dat, obs = "obs", mod = "mod", group = "model", normalise = TRUE)


# Let's replicate this...
# First, gathering the statistics...
corrcoeff_func_simp<- function(df){
  df.use<- df %>%
    drop_na(obs, mod)
  mean.obs<- mean(df.use$obs)
  mean.mod<- mean(df.use$mod)
  sd.obs<- sd(df.use$obs)
  sd.mod<- sd(df.use$mod)
  samps<- nrow(df.use)

  out<- ((1/samps)*(sum((df.use$mod - mean.mod)*(df.use$obs - mean.obs))))/(sd.obs*sd.mod)
  #out <- cor(df.use$obs, df.use$mod, use = "pairwise")
  return(out)
}

bias_func_simp<- function(df){
  df.use<- df %>%
    drop_na(obs, mod)
  out<- sd(df.use$mod)/sd(df.use$obs)
  return(out)
}

mod.dat.stats<- mod.dat %>%
  group_by(model) %>%
  nest() %>%
  mutate(., "CorrCoeff" = as.numeric(map(data, corrcoeff_func_simp)),
         "Bias" = as.numeric(map(data, bias_func_simp)))

# Now plot creation....
# Getting maxSD for plotting
maxsd<- max(mod.dat.stats$Bias, 1)

# Empty plot first
# Creating empty plot first
plot.base<- ggplot() + 
  scale_x_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5)) +
  scale_y_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5)) +
  theme_classic()

# Coeff D rays 
grad.corr.lines = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for(i in 1:length(grad.corr.lines)){
  x.vec<- c(0, maxsd*grad.corr.lines[i])
  y.vec<- c(0, maxsd*sqrt(1 - grad.corr.lines[i]^2))
  
  if(i ==1){
    coeffd.rays.df<- data.frame("Ray" = rep(1, length(x.vec)), "x" = x.vec, "y" = y.vec)
  } else {
    temp<- data.frame("Ray" = rep(i, length(x.vec)), "x" = x.vec, "y" = y.vec)
    coeffd.rays.df<- bind_rows(coeffd.rays.df, temp)
  }
}

# Add rays
plot.coeffd<- plot.base +
  geom_line(data = coeffd.rays.df, aes(x = x, y = y, group = Ray), lty = "longdash", col = "lightgray")

coeffd.labs<- coeffd.rays.df %>%
  group_by(Ray) %>%
  summarize(., 
            "x" = max(x, na.rm = TRUE), 
            "y" = max(y, na.rm = TRUE)) %>%
  data.frame()

coeffd.labs$Label<- grad.corr.lines

plot.coeffd<- plot.coeffd +
  geom_label(data = coeffd.labs, aes(x = x, y = y, label = Label), fill = "white", label.size = NA)

# SD arcs
# Need to add in SD arcs
sd.arcs<- seq(from = 0, to = maxsd, by = 0.5)

for(i in 1:length(sd.arcs)){
  x.vec<- sd.arcs[i]*cos(seq(0, pi/2, by = 0.03))
  y.vec<- sd.arcs[i]*sin(seq(0, pi/2, by = 0.03))
  
  if(i ==1){
    sd.arcs.df<- data.frame("Arc" = rep(sd.arcs[1], length(x.vec)), "x" = x.vec, "y" = y.vec)
  } else {
    temp<- data.frame("Arc" = rep(sd.arcs[i], length(x.vec)), "x" = x.vec, "y" = y.vec)
    sd.arcs.df<- bind_rows(sd.arcs.df, temp)
  }
}

# Add arcs to plot.base
plot.sd<- plot.coeffd +
  geom_line(data = sd.arcs.df, aes(x = x, y = y, group = Arc), lty = "dotted", color = "lightgray") 

# Now gamma? -- Standard deviation arcs around the reference point
gamma<- pretty(c(0, maxsd), n = 4)[-1]
gamma<- gamma[-length(gamma)]
labelpos<- seq(45, 70, length.out = length(gamma))

for(gindex in 1:length(gamma)) {
  xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + sd.r
  endcurve <- which(xcurve < 0)
  endcurve <- ifelse(length(endcurve), min(endcurve) - 1, 105)
  ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
  maxcurve <- xcurve * xcurve + ycurve * ycurve
  startcurve <- which(maxcurve > maxsd * maxsd)
  startcurve <- ifelse(length(startcurve), max(startcurve) + 1, 0)
  x.vec<- xcurve[startcurve:endcurve]
  y.vec<- ycurve[startcurve:endcurve]
  
  if(gindex ==1){
    gamma.df<- data.frame("Gamma" = rep(gamma[1], length(x.vec)), "x" = x.vec, "y" = y.vec)
  } else {
    temp<- data.frame("Gamma" = rep(gamma[gindex], length(x.vec)), "x" = x.vec, "y" = y.vec)
    gamma.df<- bind_rows(gamma.df, temp)
  }
}

gamma.df$Gamma<- factor(gamma.df$Gamma, levels = unique(gamma.df$Gamma))

# Add em
plot.gamma<- plot.sd +
  geom_line(data = gamma.df, aes(x = x, y = y, group = Gamma), lty = "solid", col = "lightgray")

# Label...
gamma.labs<- gamma.df %>%
  group_by(Gamma) %>%
  summarize("x" = mean(x, na.rm = TRUE), 
            "y" = median(y, na.rm = TRUE))

plot.gamma<- plot.gamma +
  geom_label(data = gamma.labs, aes(x = x, y = y, label = Gamma), fill = "white", label.size = NA)

# Add in reference point
plot.all<- plot.gamma +
  geom_point(aes(x = sd.r, y = 0), color = "black", size = 4)

# Add in reference points
mod.results.td<- mod.dat.stats %>%
  mutate(., "TD.X" = Bias * CorrCoeff,
         "TD.Y" = Bias * sin(acos(CorrCoeff)))

ours.td<- plot.all +
  geom_point(data = mod.results.td, aes(x = TD.X, y = TD.Y, color = model), size = 2)

TaylorDiagram(mod.dat, obs = "obs", mod = "mod", group = "model", normalise = TRUE)
taylor.diagram(ref = mod.dat$obs[mod.dat$model == "model 1"], model = mod.dat$mod[mod.dat$model == "model 1"], normalize = TRUE)
taylor.diagram(ref = mod.dat$obs[mod.dat$model == "model 2"], model = mod.dat$mod[mod.dat$model == "model 2"], normalize = TRUE, add = TRUE, col = "blue")
taylor.diagram(ref = mod.dat$obs[mod.dat$model == "model 3"], model = mod.dat$mod[mod.dat$model == "model 3"], normalize = TRUE, add = TRUE, col = "green")

