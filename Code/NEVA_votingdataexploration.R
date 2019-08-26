#####
## NEVA voting data
#####

##
library(tidyverse)
library(gtools)
##

## What do we have??
dir.dat<- read.csv("./Data/JHareDirectionalEffect.csv")
vuln.dat<- read.csv("./Data/JHareQualitativeDataResults.csv")

## Directional effect voting: 3 individuals, each with 4 votes to place into 3 distinct bins. Could we get to the individual?

## Vulnerability: Sensitivity and Exposure
head(vuln.dat)

tsens<- vuln.dat[vuln.dat$Species == "Acadian Redfish" & vuln.dat$Attribute.Category == "Sensitivity.Attribute",]
texp<- vuln.dat[vuln.dat$Species == "Acadian Redfish" & vuln.dat$Attribute.Category == "Exposure.Factor",]

# Sensitivity: 12 factors, 5 voters, 5 votes per person to place in 4 bins: 300 total votes, 25 votes per factor 
sum(tsens %>%
  summarise_at(., .vars = c("Low", "Moderate", "High", "Very.High"), sum))

sum(tsens[1,5:8])

# Exposure: 12 factors, 4 voters, 5 votes per person to place in 4 bins: 240 votes total, 20 votes per factor
sum(texp %>%
      summarise_at(., .vars = c("Low", "Moderate", "High", "Very.High"), sum))

sum(texp[1,5:8])

## Need to get from 12 factors per species for sensitivity and exposure to one overall measurement...Mean should be fine...
tsens.avg<- tsens %>%
  summarise_at(., .vars = c("Low", "Moderate", "High", "Very.High"), mean) 
texp.avg<- texp %>%
  summarise_at(., .vars = c("Low", "Moderate", "High", "Very.High"), mean)

## Reweighting to account for difference in number of voters??
tsens.wt<- 5/5
texp.wt<- 5/4
dir.wt<- 5/3

texp.avg.rwt<- round(texp.avg * texp.wt)
tsens.avg.rwt<- round(tsens.avg * tsens.wt)

## Likelihood and Prior weights
auc.range<- seq(from = 0.5, to = 1.0, length.out = 100)

# When the auc is high (near 1), we'd expect that the range of possible smooths is relatively small. So, presumably all of the possible curves are pretty good. One option is to then inflate the importance of NEVA likelihood, such that it really explores that envelope? In contrast, with a poor model, auc will be lower (0.5). The envelope around the smooth curves for the poor model will be large and we don't really want NEVA having the ability to explore those (unlikely) curves?
neva.wts<- 2*auc.range

plot(x = auc.range, y = neva.wts)
