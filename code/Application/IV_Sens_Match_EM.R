library(foreign)
library(optmatch)
library(rcbalance)
library(plyr)
library(dplyr)

rm(list=ls())

source("./rfunctions/meantabRC.R")
source("./rfunctions/build.dist.struct.nearfar.v2.R")

# Load the Data
data <- read.dta("./data/iv_match_s1.dta")
data <- data[order(data$pref_bin, decreasing = TRUE), ]
data$id <- seq(1:nrow(data))
data$treat <- data$pref_bin
n.pat <- table(data$pref_bin)
n.pat

# Variables to Match On
vars <-  c("age","comorb", "med_income",
           "income_miss", "female", "hisp", "white", "afam", "other", 
           "angus", "disability_bin", "vol_low", "vol_med", "vol_high",
           "age_gr1","age_gr2","age_gr3","age_gr4","age_gr5", "age_gr6",
           "p_cat1", "p_cat2", "p_cat3", "p_cat4", "p_cat5",
           "agenow", "yearsexp", "exp_miss", "mort_rate", "plos_rate",
           "ynel1", "ynel2","ynel3","ynel4", "ynel5","ynel6", "ynel7",
           "ynel8","ynel9", "ynel10","ynel11","ynel12","ynel13","ynel14",
           "ynel15","ynel16","ynel17", "ynel18","ynel19","ynel20","ynel21",
           "ynel22","ynel23","ynel24","ynel25","ynel26", "ynel27","ynel28",
           "ynel29","ynel30","ynel31")
far.var <- c("pref")


# Build Distance Matrix w/ Caliper to Strgth IV
t_fexact <- system.time({  
  my.dist <- build.dist.struct.nearfar(z= data$pref_bin, 
                                       X=as.matrix(data[vars]),
                                       calip.option = "propensity", 
                                       calip.cov = NULL, caliper = .08, 
                                       reverse.calip.cov = as.matrix(data[far.var]), 
                                       reverse.caliper = 1.0, 
                                       verbose=TRUE, 
                                       exact = paste(data$angus, data$hospid,
                                                     sep = '.')
  )
  
})
# time required for computation, in minutes
t_fexact[['elapsed']]/60

#define fine balance levels
l1 <- c("white")
l2 <- c(l1, "vol_low", "vol_med", "vol_high")
l3 <- c(l1, l2, "disability_bin")

# Run the Match
t_fexact <- system.time({       
  match.out <- rcbalance(my.dist, fb.list = list(l1, l2, l3),
                         treated.info = data[data$pref_bin == 1,],
                         control.info = data[data$pref_bin != 1,], 
                         exclude.treated = TRUE, tol = .001)
  
})

# time required for computation, in minutes
t_fexact[['elapsed']]/60

#Post Match Processing
t.dat <- data[data$pref_bin==1,]
c.dat <- data[data$pref_bin==0,]
match.treat <- t.dat[as.numeric(rownames(match.out$matches)),]
n <- nrow(match.treat)
# How many people dropped in match to respect within hospital contraint?
n*2
nrow(data)

# Create the Matched Data Set
match.treat$pair.id <- 1:n
match.ctrl <- c.dat[match.out$matches,]
match.ctrl$pair.id <- 1:n

match.data <- rbind(match.treat, match.ctrl)
bal_vars <-  c("id", "treat", "pref", "age","comorb", "med_income",
           "income_miss", "female", "hisp", "white", "afam", "other", 
           "angus", "disability_bin", "vol_low", "vol_med", "vol_high",
           "age_gr1","age_gr2","age_gr3","age_gr4","age_gr5", "age_gr6",
           "p_cat1", "p_cat2", "p_cat3", "p_cat4", "p_cat5",
           "agenow", "yearsexp", "exp_miss", "mort_rate", "plos_rate",
           "ynel1", "ynel2","ynel3","ynel4", "ynel5","ynel6", "ynel7",
           "ynel8","ynel9", "ynel10","ynel11","ynel12","ynel13","ynel14",
           "ynel15","ynel16","ynel17", "ynel18","ynel19","ynel20","ynel21",
           "ynel22","ynel23","ynel24","ynel25","ynel26", "ynel27","ynel28",
           "ynel29","ynel30","ynel31")

# Calculate a Balance Table
X_mat = data[bal_vars]
ids <- match.data$id
bal.tab.post <- meantabRC(X_mat, ids, exact=NULL, digits=2)
row.names(bal.tab.post) <- bal_vars[c(-1:-2)]
bal.short <- bal.tab.post[which(abs(bal.tab.post[,6])>.09),]

# Save the Results
save(bal.tab.post, bal.short, match.data, X_mat, file ="./routput/iv_match_final_em.RData")




