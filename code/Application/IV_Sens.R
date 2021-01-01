library(foreign)
library(dplyr)

rm(list=ls())

source("./rfunctions/outcome_funcs.R")
source("./rfunctions/ERSensitivityBinary.R")
source("./rfunctions/ERSensitivityHet.R")
source("./rfunctions/EffectRatio.R")
source("./rfunctions/SEPoP.R")
source("./rfunctions/ERSensitivityPop.R")
source("./rfunctions/EffModTestPop.R")

# Load the Data
load("./routput/iv_match_final_em.RData")

######## Sepsis Subgroup

## Stratify the Data
data.sep <- match.data[match.data$angus==1,]
n.p <- nrow(data.sep)/2
data.sep$pair.id <- c(seq(1: n.p), seq(1: n.p))
								
#### IV Sens Binary Outcome
sen.comp.sep <- ERSensitivityBinary(data.sep$pair.id, data.sep$complication, data.sep$pref_bin, data.sep$surg, 
                                    null=0, 
                                    Gamma = 1, 
									monotonicity = T, 
									exclusionRestriction = T, 
									effectDirection = "both", 
									Changepoint = T,
									SensitivityInterval = T, 
									nperm = 1000, 
									alternative = "two.sided", 
									continuousRelax = T)


# IV Sens Continuous Outcome
Ddose.sep = data.sep$surg[data.sep$pref_bin==1] - data.sep$surg[data.sep$pref_bin==0]									
sen.los.sep = ERSensitivity(y.los.sep, Ddose.sep, 
                             Gamma = 1,  
                             Changepoint = T, 
                             SensitivityInterval = T, 
                             nperm = 1000, 
                             alternative = "two.sided")

######## Non Sepsis Subgroup

## Stratify the Data
data.nsep <- match.data[match.data$angus==0,]
n.p <- nrow(data.nsep)/2
data.nsep$pair.id <- c(seq(1: n.p), seq(1: n.p))

#### IV Sens Binary Outcome
sen.comp.nsep <- ERSensitivityBinary(data.nsep$pair.id, data.nsep$complication, data.nsep$pref_bin, data.nsep$surg, 
                                    null=0, 
                                    Gamma = 1, 
									monotonicity = T, 
									exclusionRestriction = T, 
									effectDirection = "both", 
									Changepoint = T,
									SensitivityInterval = T, 
									nperm = 1000, 
									alternative = "two.sided", 
									continuousRelax = T)

# IV Sens Continuous Outcome
Ddose.nsep = data.nsep$surg[data.nsep$pref_bin==1] - data.nsep$surg[data.nsep$pref_bin==0]									
sen.los.nsep = ERSensitivity(y.los.nsep, Ddose.nsep, 
                             Gamma = 1,  
                             Changepoint = T, 
                             SensitivityInterval = T, 
                             nperm = 1000, 
                             alternative = "two.sided")

# Save the Results
save(sen.comp.sep,sen.los.sep,     
     sen.comp.nsep,sen.los.nsep,
     file="./routput/sens_results.RData")
     
    