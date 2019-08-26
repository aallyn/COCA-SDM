# Experimenting with SDM Bayes
library(mgcv)
library(tidyverse)
library(gridExtra)

func.path<- "./Code/"
source(paste(func.path, "sdm_bayesianupdating_func.R", sep = ""))

#############
### Andy's examples
x.vec<- seq(from = -5, to = 5, length.out = 500)
base<- dnorm(x.vec, mean = 0, sd = 1)
fut<- dnorm(x.vec, mean = 1, sd = 1)

fut.means<- seq(-3, 3, length.out = 20)
means<- rep(NA, length(fut.means))
sds<- rep(NA, length(fut.means))

for(i in 1:length(fut.means)){
  out<- SDMbayes(x.vec, 0, 1, fut.means[i], 1, 0, 3, 9, 0, 0, 3, 9)
  means[i]<- x.vec%*%out/sum(out)
  #sds[i]<- sd(out)
  print(fut.means[i])
}

plot(fut.means, means) 

sdm.bayes1<- SDMbayes(x.vec, 0, 1, 1, 1, 0, 3, 9, 0, 0, 3, 9)
sdm.bayes2<- SDMbayes(x.vec, 0, 1, 1, 1, 0, 3, 9, 0, 2, 5, 5)
sdm.bayes3<- SDMbayes(x.vec, 0, 1, 1, 1, 0, 3, 9, 0, 4, 4, 4)
sdm.bayes4<- SDMbayes(x.vec, 0, 1, 1, 1, 9, 3, 0, 4, 4, 4, 0)

plot.dat<- data.frame("Sample" = c(rep("SDM.Base", length(base)), rep("SDM.Future", length(fut)), rep("Pos.Vh", length(sdm.bayes1)), rep("Pos.H", length(sdm.bayes2)), rep("Pos.M", length(sdm.bayes3)), rep("Neg.M", length(sdm.bayes4))), "X" = rep(x.vec, 6), "Value" = c(base, fut, sdm.bayes1, sdm.bayes2, sdm.bayes3, sdm.bayes4))
plot.dat$Sample<- factor(plot.dat$Sample, levels = c("SDM.Base", "SDM.Future", "Pos.Vh", "Pos.H", "Pos.M", "Neg.M"))

out.plot<- ggplot(plot.dat, aes(x = X, y = Value, group = Sample)) + 
  geom_line(aes(color = Sample), alpha = 0.75) +
  scale_fill_manual(name = "Sample", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33'), labels = c("SDM.Base", "SDM.Future", "Pos.Vh", "Pos.H", "Pos.M", "Neg.M")) +
  theme_bw()
out.plot

## This seems promising and is the same thing we are getting in Matlab with Andy's original code, so the translation seems to have worked. 

###############
### Function behavior...assessing how the function behaves at limits (i.e., going from a bmn of 0 to fmn of 1, or bmn of 1 to fmn of 0) across different voting situations. Note:: This change (0-1 and 1-0) is interesting on the RESPONSE scale, but we are actually working on the link scale...to get something equivalent, we could use -5 and 3? Why?
logit_func<- function(x) {
  exp(x)/(1+exp(x))
}

logit_func(-5)
logit_func(3)

## Starting simple, lets loop over a vector of differences between those values
# baseline mean possibilities
bmn<- seq(from = -5, to = 3, by = 0.5)

# future mean possibilities
fmn<- seq(from = 3, to = -5, by = -0.5)

base.fut.poss<- data.frame("Scenario.Base.Fut" = c(rep("Increasing P(Presence)", 8), "No Change", rep("Decreasing P(Presence)", 8)), "Base.Mean" = bmn, "Fut.Mean" = fmn, "Presence.Change" = fmn-bmn)

# What did that just do? -- creates a dataframe ranging in presence changes (base to future) from 8 to -8
base.fut.poss

# Now the loop -- this would mean keeping the directional effect voting and the vulnerability voting constant.
SDMBayes.loop<- vector("list", nrow(base.fut.poss))
nevaD.neg1<- 0
nevaD.neut1<- 0.1
nevaD.pos1<- 0
nevaV.low1<- 0.1
nevaV.mod1<- 0
nevaV.high1<- 0
nevaV.vhigh1<- 0
x.vec1<- seq(from = -10, to = 10, length.out = 500)

for(i in 1:nrow(base.fut.poss)){
  SDMBayes.loop[[i]]<- SDMbayes(x.vec1, 0, 1, base.fut.poss$Fut.Mean[i], 1, nevaD.neg1, nevaD.neut1, nevaD.pos1, nevaV.low1, nevaV.mod1, nevaV.high1, nevaV.vhigh1)
  names(SDMBayes.loop)[i]<- paste("Change_", base.fut.poss$Presence.Change[i], sep = "")
}

# Mean from each?
test_fun<- function(x, x.vec.use = x.vec){
  out<- x%*%x.vec.use/sum(x)
  return(out)
}
SDMBayes.loop.means<- data.frame("Presence.Change" = base.fut.poss$Presence.Change, "Mean" = unlist(lapply(SDMBayes.loop, test_fun)))
SDMBayes.loop.means

plot.dat<- data.frame("Sample" = c(rep("SDM.1v1", 17), rep("SDM.NEVA", 17)), "Xaxis.Link" = rep(base.fut.poss$Fut.Mean, 2), "Xaxis.Resp" = rep(base.fut.poss$Fut.Mean, 2),  "LinkValue" = c(base.fut.poss$Fut.Mean, SDMBayes.loop.means$Mean), "RespValue" = c(ilink(base.fut.poss$Fut.Mean), ilink(SDMBayes.loop.means$Mean)))
plot.dat$Sample<- factor(plot.dat$Sample, levels = c("SDM.1v1", "SDM.NEVA"))

link.plot<- ggplot(plot.dat, aes(x = Xaxis.Link, y = LinkValue, group = Sample)) +
  geom_line(aes(color = Sample), alpha = 0.75) +
  scale_fill_manual(name = "Sample", values = c('#e41a1c','#377eb8'), labels = c("SDM.1v1", "SDM.NEVA")) +
  ggtitle("LinkScale") +
  theme_bw()
  
resp.plot<- ggplot(plot.dat, aes(x = Xaxis.Resp, y = RespValue, group = Sample)) +
  geom_line(aes(color = Sample), alpha = 0.75) +
  scale_fill_manual(name = "Sample", values = c('#e41a1c','#377eb8'), labels = c("SDM.1v1", "SDM.NEVA")) +
  ggtitle("ResponseScale") +
  theme_bw()

grid.arrange(link.plot, resp.plot, ncol = 2)

















# Now, create some possibilities for the Directional effect votes and the Vulnerability votes
# nevaD has 12 possible votes...
nevaD.poss<- data.frame("Scenario.Dir" = c("Neg.Vh", "Neg.M", "Neg.L", "Neut.Vh", "Neut.M", "Neut.L", "Pos.Vh", "Pos.M", "Pos.L"), "Negative" = c(12, 9, 6, 0, 2, 3, 0, 0, 3), "Neutral" = c(0, 3, 3, 12, 8, 6, 0, 3, 3), "Positive" = c(0, 0, 3, 0, 2, 3, 12, 9, 6))
nevaD.poss

#nevaV (considering 24 possible votes)
nevaV.poss<- data.frame("Scenario.Vuln" = c("Low.Vh", "Low.M", "Low.L", "Mod.Vh", "Mod.M", "Mod.L", "High.Vh", "High.M", "High.L", "VHigh.Vh", "VHigh.M", "VHigh.L"), "Low" = c(24, 20, 10, 0, 3, 6, 0, 0, 2, 0, 0, 2), "Mod" = c(0, 4, 6, 24, 18, 10, 0, 3, 6, 0, 0, 6), "High" = c(0, 0, 6, 0, 3, 6, 24, 18, 10, 0, 4, 6) , "VHigh" = c(0, 0, 2, 0, 0, 2, 0, 3, 6, 24, 20, 10))
nevaV.poss

# Expand these three things: presence change, directional effect voting, vulnerability voting -- ALL combinations
scen.combo<- expand.grid("Scenario.Dir" = nevaD.poss$Scenario.Dir, "Scenario.Vuln" = nevaV.poss$Scenario.Vuln, "Presence.Change" = base.fut.poss$Presence.Change)
scen.combo<- scen.combo %>%
  left_join(., nevaD.poss, by = "Scenario.Dir") %>%
  left_join(., nevaV.poss, by = "Scenario.Vuln") %>%
  left_join(., base.fut.poss, by = "Presence.Change")
View(scen.combo)

# Map SDMBayes function to each voting line
x.vec<- seq(from = -10, to = 10, length.out = 500)

scen.combo<- scen.combo %>%
  mutate(., "SDMBayes" = pmap(list(x = list(x.vec), bmn = Base.Mean, bsd = list(1), fmn = Fut.Mean, fsd = list(1), nevaD.neg = Negative, nevaD.neut = Neutral, nevaD.pos = Positive, nevaV.low = Low, nevaV.mod = Mod, nevaV.high = High, nevaV.vhigh = VHigh), SDMbayes))

# What the heck did that do?? - Ran our SDMbayes function on each line of scen.combo, as each line has all the necessary arguments to run the SDMbayes function (think of it as an input file).

# Just a quick check -- pull out the third row. Use those inputs and run SDMbayes outside of map. Plot results together.
check.input<- data.frame(scen.combo[3,-14])
x.vec<- seq(from = -10, to = 10, length.out = 500)
base<- dnorm(x.vec, mean = check.input$Base.Mean, sd = 1)
fut<- dnorm(x.vec, mean = check.input$Fut.Mean, sd = 1)

sdm.bayes.man<- SDMbayes(x.vec, check.input$Base.Mean, 1, check.input$Fut.Mean, 1, check.input$Negative, check.input$Neutral, check.input$Positive, check.input$Low, check.input$Mod, check.input$High, check.input$VHigh)
sdm.bayes.map<- scen.combo$SDMBayes[[3]]

plot.dat<- data.frame("Sample" = c(rep("SDM.Base", length(base)), rep("SDM.Future", length(fut)), rep("Manual", length(sdm.bayes.man)), rep("Mapped", length(sdm.bayes.map))), "X" = rep(x.vec, 4), "Value" = c(base, fut, sdm.bayes.man, sdm.bayes.map))
plot.dat$Sample<- factor(plot.dat$Sample, levels = c("SDM.Base", "SDM.Future", "Manual", "Mapped"))

out.plot<- ggplot(plot.dat, aes(x = X, y = Value, group = Sample)) + 
  geom_line(aes(color = Sample), alpha = 0.75) +
  scale_fill_manual(name = "Sample", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3'), labels = c("SDM.Base", "SDM.Future", "Manual", "Mapped")) +
  theme_bw()
out.plot

## So, the mapping seems to be working. The plot is a bit weird --- maybe because of the link vs. response scale issue?

# Now, for the plot Andy was envisioning...I think we want to pull out the no change scenarios...and then subtract each of the SDMBayes fits from this no change scenario. This would basically give us the influence of the directional effect/vulnerability rank?
scen.combo.base<- scen.combo %>%
  dplyr::filter(., Presence.Change == 0.0) %>%
  mutate(., "Merge.Col" = paste(Scenario.Dir, "_", Scenario.Vuln, sep = "")) %>%
  dplyr::select(., SDMBayes, Merge.Col)
names(scen.combo.base)[1]<- "SDMBayes.Base"

# Now, let's add that new "SDMBayes.Base" column back in and this will allow us to subtract the Base from each scenario?
scen.combo<- scen.combo %>%
  dplyr::filter(., Presence.Change != 0.0) %>%
  mutate(., "Merge.Col" = paste(Scenario.Dir, "_", Scenario.Vuln, sep = "")) %>%
  left_join(., scen.combo.base, by = "Merge.Col") %>%
  as_tibble()

# Write a difference function to map to each row
diff_func<- function(New, Base){
  out<- New-Base
  return(out)
}

# Apply it
scen.combo<- scen.combo %>%
  mutate(., "SDMBayes.Diff" = map2(SDMBayes, SDMBayes.Base, diff_func))

## Visualizing...
# Select a directional effect and vulnerability scenario to visualize. Note -- I actually don't think this is what we want, really we want difference in the pdfs??
temp<- scen.combo %>%
  dplyr::filter(Scenario.Dir == "Pos.Vh" & Scenario.Vuln == "Low.Vh") %>%
  mutate(., mean.diff = map(SDMBayes.Diff, mean),
         sd = map(SDMBayes.Diff, sd)) %>%
  dplyr::select(., -SDMBayes, -SDMBayes.Base, -SDMBayes.Diff) %>%
  unnest() %>%
  data.frame()

p<- ggplot(temp, aes(x=Presence.Change, y=mean.diff)) + 
  geom_point(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.diff-sd, ymax=mean.diff+sd), width=.2,
                position=position_dodge(.9)) 
p

# Nope...not good...difference in pdfs? mean(a-b) = meana - meanb, var(a-b) = var(a) + var(b)
norm_diff_func<- function(New, Base, Stat){
  if(Stat == "mean"){
    out<- mean(New) - mean(Base)
    return(out)
  }
  
  if(Stat == "sd"){
    out<- sqrt(var(New) + var(Base))
    return(out)
  }
}

# Try that
temp<- scen.combo %>%
  dplyr::filter(Scenario.Dir == "Pos.Vh" & Scenario.Vuln == "Low.Vh") %>%
  mutate(., mean.diff.pdf = pmap(list(New = SDMBayes, Base = SDMBayes.Base, Stat = "mean"), norm_diff_func),
         sd.pdf = pmap(list(New = SDMBayes, Base = SDMBayes.Base, Stat = "sd"), norm_diff_func)) %>%
  dplyr::select(., -SDMBayes, -SDMBayes.Base, -SDMBayes.Diff) %>%
  unnest() %>%
  data.frame()

p<- ggplot(temp, aes(x=Presence.Change, y=mean.diff.pdf)) + 
  geom_point(stat="identity", color="black", 
             position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.diff.pdf-sd.pdf, ymax=mean.diff.pdf+sd.pdf), width=.2,
                position=position_dodge(.9)) 
p

########
## What about with real observations?
# Load in our model fits to get bmn and bsd
all.dat<- readRDS("./Data/sdm.projections.SST_01172018.rds") 

# Pull out a row
temp.dat<- all.dat[1,]

ilink <- family(all.dat$Model.Fitted[[1]])$linkinv

# Draw n.samps from normal distribution with mean  = pred.dat mu 
# and sd = pred.dat se. 
# baseline
base.vec.link<- rnorm(500, mean = temp.dat$Projections[[1]]$Baseline[1], sd = temp.dat$Projections.p.se[[1]]$Baseline[1])
base.vec<- ilink(base.vec.link)

# future
fut.vec.link<- rnorm(500, mean = temp.dat$Projections[[1]]$`2055`[1], sd = temp.dat$Projections.p.se[[1]]$`2055`[1])
fut.vec<- ilink(fut.vec.link)

# Propose x.vector values (link scale, +/- 5SD from the mean)
x.vec<- seq(from = temp.dat$Projections[[1]]$Baseline[1] + 5*temp.dat$Projections.p.se[[1]]$Baseline[1], to = temp.dat$Projections[[1]]$Baseline[1] - 5*temp.dat$Projections.p.se[[1]]$Baseline[1], length.out = 500)
base.r<- ilink(dnorm(x.vec, mean = temp.dat$Projections[[1]]$Baseline[1], sd = temp.dat$Projections.p.se[[1]]$Baseline[1]))
fut.r<- ilink(dnorm(x.vec, mean = temp.dat$Projections[[1]]$`2055`[1], sd = temp.dat$Projections.p.se[[1]]$`2055`[1]))

bmn = temp.dat$Projections[[1]]$Baseline[1]
bsd = temp.dat$Projections.p.se[[1]]$Baseline[1]
fmn = temp.dat$Projections[[1]]$`2055`[1]
fsd = temp.dat$Projections.p.se[[1]]$`2055`[1]
nevaD = c(9, 3, 0)
nevaV = c(90, 10, 0, 0)

sdm.bayes1.r<- ilink(SDMbayes(x.vec, bmn, bsd, fmn, fsd, 0, 3, 9, 0, 1, 10, 20))
sdm.bayes2.r<- ilink(SDMbayes(x.vec, bmn, bsd, fmn, fsd, 0, 3, 9, 0, 1, 20, 10))
sdm.bayes3.r<- ilink(SDMbayes(x.vec, bmn, bsd, fmn, fsd, 0, 3, 9, 0, 6, 15, 10))
sdm.bayes4.r<- ilink(SDMbayes(x.vec, bmn, bsd, fmn, fsd, 9, 3, 0, 0, 6, 15, 10))

plot.dat<- data.frame("Sample" = c(rep("SDM.Base", length(base.r)), rep("SDM.Future", length(fut.r)), rep("Pos.Vh", length(sdm.bayes1.r)), rep("Pos.H", length(sdm.bayes2.r)), rep("Pos.M", length(sdm.bayes3.r)), rep("Neg.M", length(sdm.bayes4.r))), "X" = rep(ilink(x.vec), 6), "Value" = c(base.r, fut.r, sdm.bayes1.r, sdm.bayes2.r, sdm.bayes3.r, sdm.bayes4.r))
plot.dat$Sample<- factor(plot.dat$Sample, levels = c("SDM.Base", "SDM.Future", "Pos.Vh", "Pos.H", "Pos.M", "Neg.M"))

out.plot<- ggplot(plot.dat, aes(x = X, y = Value, group = Sample)) + 
  ylim(c(0,1)) +
  geom_line(aes(color = Sample), alpha = 0.75) +
  scale_fill_manual(name = "Sample", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33'), labels = c("SDM.Base", "SDM.Future", "Pos.Vh", "Pos.H", "Pos.M", "Neg.M")) +
  theme_bw()
out.plot

# NICE! Why are we at a minimum of 0.5 on the y axis??

# Now, how do things work when we have a spatial gradient in species increasing and decreasing across different votes. Imagine two pixels, one shifted up and one shifted down, and then how these respond to different voting ---

# Good example of species P(presence) increasing or decreasing? Longfin squid southern NS shelf increase
temp.dat<- all.dat[75,]

diff<- ilink(temp.dat$Projections[[1]]$`2055`) - ilink(temp.dat$Projections[[1]]$`Baseline`)
diff.max.ind<- which.max(diff) # Row 87
diff.max<- diff[diff.max.ind]

bmn.inc = temp.dat$Projections[[1]]$Baseline[diff.max.ind]
bsd.inc = temp.dat$Projections.p.se[[1]]$Baseline[diff.max.ind]
fmn.inc = temp.dat$Projections[[1]]$`2055`[diff.max.ind]
fsd.inc = temp.dat$Projections.p.se[[1]]$`2055`[diff.max.ind]

# nevaD has 12 possible votes...
nevaD.poss<- data.frame("Scenario.Dir" = c("Neg.Vh", "Neg.M", "Neg.L", "Neut.Vh", "Neut.M", "Neut.L", "Pos.Vh", "Pos.M", "Pos.L"), "Negative" = c(12, 9, 6, 0, 2, 3, 0, 0, 3), "Neutral" = c(0, 3, 3, 12, 8, 6, 0, 3, 3), "Positive" = c(0, 0, 3, 0, 2, 3, 12, 9, 6))

#nevaV has 230 possible votes
nevaV.poss<- data.frame("Scenario.Vuln" = c("Low.Vh", "Low.M", "Low.L", "Mod.Vh", "Mod.M", "Mod.L", "High.Vh", "High.M", "High.L", "VHigh.Vh", "VHigh.M", "VHigh.L"), "Low" = c(24, 20, 10, 0, 3, 6, 0, 0, 2, 0, 0, 2), "Mod" = c(0, 4, 6, 24, 18, 10, 0, 3, 6, 0, 0, 6), "High" = c(0, 0, 6, 0, 3, 6, 24, 18, 10, 0, 4, 6) , "VHigh" = c(0, 0, 2, 0, 0, 2, 0, 3, 6, 24, 20, 10))


# Table of all possible combinations...
scen.combo<- expand.grid("Scenario.Dir" = nevaD.poss$Scenario.Dir, "Scenario.Vuln" = nevaV.poss$Scenario.Vuln)
scen.combo<- scen.combo %>%
  left_join(., nevaD.poss, by = "Scenario.Dir") %>%
  left_join(., nevaV.poss, by = "Scenario.Vuln") %>%
  as_tibble()

# Map SDMBayes function to each line...a cell going from absent to present
x.vec<- seq(from = -6, to = 4, length.out = 500)
bmn.inc<- -5
bsd.inc<- 0.2
fmn.inc<- 3
fsd.inc<- 0.2

scen.combo<- scen.combo %>%
  mutate(., "SDMBayes.Inc" = pmap(list(x = list(x.vec), bmn = list(bmn.inc), bsd = list(bsd.inc), fmn = list(fmn.inc), fsd = list(fsd.inc), nevaD.neg = Negative, nevaD.neut = Neutral, nevaD.pos = Positive, nevaV.low = Low, nevaV.mod = Mod, nevaV.high = High, nevaV.vhigh = VHigh), SDMbayes))


# Findings --- MOSTLY NAs


# Map SDMBayes function to each line...a cell with no change (1:1)
x.vec<- seq(from = -6, to = 4, length.out = 500)
bmn.neut<- 3
bsd.neut<- 0.2
fmn.neut<- 3
fsd.neut<- 0.2

scen.combo<- scen.combo %>%
  mutate(., "SDMBayes.Neut" = pmap(list(x = list(x.vec), bmn = list(bmn.neut), bsd = list(bsd.neut), fmn = list(fmn.neut), fsd = list(fsd.neut), nevaD.neg = Negative, nevaD.neut = Neutral, nevaD.pos = Positive, nevaV.low = Low, nevaV.mod = Mod, nevaV.high = High, nevaV.vhigh = VHigh), SDMbayes))


# All fine?
plot(ilink(x.vec), ilink(scen.combo$SDMBayes.Neut[[1]])) 
plot(ilink(x.vec), ilink(scen.combo$SDMBayes.Neut[[88]])) 


# Map SDMBayes function to each line...a cell going from present to absent
x.vec<- seq(from = -6, to = 4, length.out = 500)
bmn.neg<- 3
bsd.neg<- 0.2
fmn.neg<- -5
fsd.neg<- 0.2

scen.combo<- scen.combo %>%
  mutate(., "SDMBayes.Neg" = pmap(list(x = list(x.vec), bmn = list(bmn.neg), bsd = list(bsd.neg), fmn = list(fmn.neg), fsd = list(fsd.neg), nevaD.neg = Negative, nevaD.neut = Neutral, nevaD.pos = Positive, nevaV.low = Low, nevaV.mod = Mod, nevaV.high = High, nevaV.vhigh = VHigh), SDMbayes))


# Findings --- cannot have a large increase with a negative species (vh certainty) and low vuln (vh certainty)
plot(ilink(x.vec), ilink(scen.combo$SDMBayes.Neg[[1]])) 
plot(ilink(x.vec), ilink(scen.combo$SDMBayes.Neg[[88]])) 
