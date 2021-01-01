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

## Stratify Data
data.sep <- match.data[match.data$angus==1,]
n.p <- nrow(data.sep)/2
data.sep$pair.id <- c(seq(1: n.p), seq(1: n.p))

# Take Outcome Matched Pair Differences
y.comp.sep <- data.sep$complication[data.sep$pref_bin==1] - data.sep$complication[data.sep$pref_bin==0] 
y.los.sep <- data.sep$los[data.sep$pref_bin==1] - data.sep$los[data.sep$pref_bin==0]
Ddose.sep = data.sep$surg[data.sep$pref_bin==1] - data.sep$surg[data.sep$pref_bin==0]	

vars <-  c("age","comorb", "med_income",
           "income_miss", "female", "hisp", "white", "afam", 
           "disability_bin", "vol_low", "vol_med", 
           "age_gr1","age_gr2","age_gr3","age_gr4","age_gr5",
           "p_cat1", "p_cat2", "p_cat3", "p_cat4", 
           "agenow", "yearsexp", "exp_miss", "mort_rate", "plos_rate",
           "ynel1", "ynel2","ynel3","ynel4", "ynel5","ynel6", "ynel7",
           "ynel8","ynel9", "ynel10","ynel11","ynel12","ynel13","ynel14",
           "ynel15","ynel16","ynel17", "ynel18","ynel19","ynel20","ynel21",
           "ynel22","ynel23","ynel24","ynel25","ynel26", "ynel27","ynel28",
           "ynel29","ynel30", "pair.id")
           
# Linear Model for SE
X.mat.sep <- data.sep[vars]
X.avg.sep <- X.mat.sep %>% group_by(pair.id) %>% summarise_all(funs(mean))
X.avg.sep <- as.data.frame(X.avg.sep[,-1])

# Sens Analysis
sen.comp.sep.lin.het = ERSensitivityHet(y.comp.sep, Ddose.sep, X = as.matrix(X.avg.sep), 
                                        Gamma = 1,  
                                        Changepoint = T, 
                                        SensitivityInterval = T, 
                                        nperm = 1000, 
                                        alternative = "two.sided")
                                        
sen.los.sep.lin.het = ERSensitivityHet(y.los.sep, Ddose.sep, X = as.matrix(X.avg.sep), 
                                       Gamma = 1,  
                                       Changepoint = T, 
                                       SensitivityInterval = T, 
                                       nperm = 1000, 
                                       alternative = "two.sided")

# Matching Based Sens Analysis
# Pair to Pair Match
Xt.sep <- data.sep[data.sep$pref_bin==1, vars]
Xc.sep <- data.sep[data.sep$pref_bin==0, vars]

# Run an NBP Match
ind.sep = indexPoP(Xt.sep, Xc.sep)

# Sens Analysis
sen.comp.sep.het <-  ERSensitivityPoP(y.comp.sep, Ddose.sep, index = ind.sep, 
                                     Gamma = 1,  
                                     Changepoint = T, 
                                     SensitivityInterval = T, 
                                     nperm = 5000, 
                                     alternative = "two.sided")
                                     
sen.los.sep.het <-  ERSensitivityPoP(y.los.sep, Ddose.sep, index = ind.sep, 
                                     Gamma = 1,  
                                     Changepoint = T, 
                                     SensitivityInterval = T, 
                                     nperm = 5000, 
                                     alternative = "two.sided")

# Calc p-values
sen.comp.sep.het.pval <-  EffModTestPoP(y.comp.sep, Ddose.sep, index = ind.sep, 
                                     nperm = 5000)
                                     
sen.los.sep.het.pval <-  EffModTestPoP(y.los.sep, Ddose.sep, index = ind.sep, 
                                     nperm = 5000)

######## Non Sepsis Subgroup

# Stratify the Daata
data.nsep <- match.data[match.data$angus==0,]
n.p <- nrow(data.nsep)/2
data.nsep$pair.id <- c(seq(1: n.p), seq(1: n.p))

# Take Outcome Matched Pair Differences
y.comp.nsep <- data.nsep$complication[data.nsep$pref_bin==1] - data.nsep$complication[data.nsep$pref_bin==0] 
y.los.nsep <- data.nsep$los[data.nsep$pref_bin==1] - data.nsep$los[data.nsep$pref_bin==0]
Ddose.nsep = data.nsep$surg[data.nsep$pref_bin==1] - data.nsep$surg[data.nsep$pref_bin==0]

# Linear Model for SE
X.mat.nsep <- data.nsep[vars]
X.avg.nsep <- X.mat.nsep %>% group_by(pair.id) %>% summarise_all(funs(mean))
X.avg.nsep <- as.data.frame(X.avg.nsep[,-1])

# Sens Analysis
sen.comp.nsep.lin.het = ERSensitivityHet(y.comp.nsep, Ddose.nsep, X = as.matrix(X.avg.nsep), 
                                         Gamma = 1,  
                                         Changepoint = T, 
                                         SensitivityInterval = T, 
                                         nperm = 1000, 
                                         alternative = "two.sided")
                                         
sen.los.nsep.lin.het = ERSensitivityHet(y.los.nsep, Ddose.nsep, X = as.matrix(X.avg.nsep), 
                                        Gamma = 1,  
                                        Changepoint = T, 
                                        SensitivityInterval = T, 
                                        nperm = 1000, 
                                        alternative = "two.sided")

# Matching Based Sens Analysis
# Pair to Pair Method
Xt.nsep <- data.nsep[data.nsep$pref_bin==1, vars]
Xc.nsep <- data.nsep[data.nsep$pref_bin==0, vars]

# Run an NBP Match
ind.nsep = indexPoP(Xt.nsep, Xc.nsep)

#Sens Analysis
sen.comp.nsep.het <-  ERSensitivityPoP(y.comp.nsep, Ddose.nsep, index = ind.nsep, 
                                     Gamma = 1,  
                                     Changepoint = T, 
                                     SensitivityInterval = T, 
                                     nperm = 5000, 
                                     alternative = "two.sided")
                                     
sen.los.nsep.het <-  ERSensitivityPoP(y.los.nsep, Ddose.nsep, index = ind.nsep, 
                                     Gamma = 1,  
                                     Changepoint = T, 
                                     SensitivityInterval = T, 
                                     nperm = 5000, 
                                     alternative = "two.sided")

# Calc p-values
sen.comp.nsep.het.pval <-  EffModTestPoP(y.comp.nsep, Ddose.nsep, index = ind.nsep, 
                                     nperm = 5000)                                    
sen.los.nsep.het.pval <-  EffModTestPoP(y.los.nsep, Ddose.nsep, index = ind.nsep, 
                                     nperm = 5000)

sen.los.sep.lin.het <- sen.los.nep.lin.het

# Save the Results
save(sen.comp.sep.lin.het, sen.los.sep.lin.het,
     sen.comp.sep.het, sen.los.sep.het,
     sen.comp.sep.het.pval, sen.los.sep.het.pval,
     sen.comp.nsep.lin.het, sen.los.nsep.lin.het, 
     sen.comp.nsep.het, sen.los.nsep.het,
     sen.comp.nsep.het.pval, sen.los.nsep.het.pval,
     file="./routput/sens_results_effmod.RData")
     
    