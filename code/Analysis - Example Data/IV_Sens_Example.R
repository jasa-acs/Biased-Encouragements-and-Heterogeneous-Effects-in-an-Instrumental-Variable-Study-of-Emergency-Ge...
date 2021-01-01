library(dplyr)


rm(list=ls())
data <- read.csv("./data/ExampleData.csv")

#### IV Estimation Via Effect Ratio
source("./rfunctions/outcome_funcs.R")
iv.out <- gen_eff_ratio(data$Y, data$D, data$Z, 
				          data$Pair, alphaLevel = 0.05, null=0)
iv.out

## IV Sensitivity Analysis
source("./rfunctions/ERSensitivityBinary.R")
source("./rfunctions/ERSensitivityHet.R")
source("./rfunctions/SEPoP.R")
source("./rfunctions/ERSensitivityPop.R")
source("./rfunctions/EffModTestPop.R")
source("./rfunctions/ERSensitivity.R")

y.dif <- data$Y[data$Z==1] - data$Y[data$Z==0]
Ddose = data$D[data$Z==1] - data$D[data$Z==0]	

# Sensitivty Analysis w/o Effect Modification 									
sen.sep = ERSensitivity(y.dif, Ddose, 
                             Gamma = 1,  
                             Changepoint = T, 
                             SensitivityInterval = T, 
                             nperm = 1000, 
                             alternative = "two.sided")
                             
                             									
## IV Sensitivity Analysis Under Effect Modification
## Linear SE
vars <-  c("X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10", "Pair")         
X.mat <- data[vars]
X.avg <- X.mat %>% group_by(Pair) %>% summarise_all(list(mean))
X.avg <- as.data.frame(X.avg.[,-1])

# Run Sens Analysis
sen.lin.het = ERSensitivityHet(y.dif, Ddose, X = as.matrix(X.avg), 
								Gamma = 1,  
								Changepoint = T, 
								SensitivityInterval = T, 
								nperm = 1000, 
								alternative = "two.sided")
									
									
## Pair to Pair Match SE
Xt <- data[data$Z==1, vars]
Xc <- data[data$Z==0, vars]

# Run an NBP Match
ind.sep = indexPoP(Xt, Xc)

# Run Sens Analysis
sen.p2p.het <-  ERSensitivityPoP(y.dif, Ddose, index = ind.sep, 
                                     Gamma = 1,  
                                     Changepoint = T, 
                                     SensitivityInterval = T, 
                                     nperm = 5000, 
                                     alternative = "two.sided")
                                     
# Omnibus p-value Test
sen.het.pval <-  EffModTestPoP(y.dif, Ddose, index = ind.sep,
                                     nperm = 5000)
