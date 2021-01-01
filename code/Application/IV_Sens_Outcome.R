
rm(list=ls())

# Load the Effect Ratio Estimation Code
source("./outcome_funcs.R")

# Load the Matched Data
load("./routput/iv_match_final_em.RData")


######## Sepsis Subgroup
data.sep <- match.data[match.data$angus==1,]
n.p <- nrow(data.sep)/2
data.sep$pair.id <- c(seq(1: n.p), seq(1: n.p))

#### ITT Estimates
itt.died.sep <- lm(died ~ pref_bin, data=data.sep)
t1.died.sep <- c(itt.died.sep$coef[2], confint(itt.died.sep)[2,1], confint(itt.died.sep)[2,2], summary(itt.died.sep)$coef[2,4])

itt.comp.sep <- lm(complication ~ pref_bin, data=data.sep)
t1.comp.sep <- c(itt.comp.sep$coef[2], confint(itt.comp.sep)[2,1], confint(itt.comp.sep)[2,2], summary(itt.comp.sep)$coef[2,4])

itt.los.sep <- lm(los ~ pref_bin, data=data.sep)
t1.los.sep <- c(itt.los.sep$coef[2], confint(itt.los.sep)[2,1], confint(itt.los.sep)[2,2], summary(itt.los.sep)$coef[2,4])

#### IV Estimates
gef.died.sep <- gen_eff_ratio(data.sep$died, data.sep$surg, data.sep$pref, 
				          data.sep$pair.id, alphaLevel = 0.05, null=0)
gef.died.sep <- c(gef.died.sep$pe, gef.died.sep$ci.l, gef.died.sep$ci.u)

gef.comp.sep <- gen_eff_ratio(data.sep$complication, data.sep$surg, data.sep$pref, 
				          data.sep$pair.id, alphaLevel = 0.05, null=0)
gef.comp.sep <- c(gef.comp.sep$pe, gef.comp.sep$ci.l, gef.comp.sep$ci.u)

gef.los.sep <- gen_eff_ratio(data.sep$los, data.sep$surg, data.sep$pref, 
				          data.sep$pair.id, alphaLevel = 0.05, null=0)
gef.los.sep <- c(gef.los.sep$pe, gef.los.sep$ci.l, gef.los.sep$ci.u)


######## Non Sepsis Subgroup
data.nsep <- match.data[match.data$angus==0,]
n.p <- nrow(data.nsep)/2
data.nsep$pair.id <- c(seq(1: n.p), seq(1: n.p))

#### ITT Estimates
itt.died.nsep <- lm(died ~ pref_bin, data=data.nsep)
t1.died.nsep <- c(itt.died.nsep$coef[2], confint(itt.died.nsep)[2,1], confint(itt.died.nsep)[2,2], summary(itt.died.nsep)$coef[2,4])

itt.comp.nsep <- lm(complication ~ pref_bin, data=data.nsep)
t1.comp.nsep <- c(itt.comp.nsep$coef[2], confint(itt.comp.nsep)[2,1], confint(itt.comp.nsep)[2,2], summary(itt.comp.nsep)$coef[2,4])

itt.los.nsep <- lm(los ~ pref_bin, data=data.nsep)
t1.los.nsep <- c(itt.los.nsep$coef[2], confint(itt.los.nsep)[2,1], confint(itt.los.nsep)[2,2], summary(itt.los.nsep)$coef[2,4])

#### IV Estimates
gef.died.nsep <- gen_eff_ratio(data.nsep$died, data.nsep$surg, data.nsep$pref, 
				          data.nsep$pair.id, alphaLevel = 0.05, null=0)
gef.died.nsep <- c(gef.died.nsep$pe, gef.died.nsep$ci.l, gef.died.nsep$ci.u)

gef.comp.nsep <- gen_eff_ratio(data.nsep$complication, data.nsep$surg, data.nsep$pref, 
				          data.nsep$pair.id, alphaLevel = 0.05, null=0)
gef.comp.nsep <- c(gef.comp.nsep$pe, gef.comp.nsep$ci.l, gef.comp.nsep$ci.u)

gef.los.nsep <- gen_eff_ratio(data.nsep$los, data.nsep$surg, data.nsep$pref, 
				          data.nsep$pair.id, alphaLevel = 0.05, null=0)
gef.los.nsep <- c(gef.los.nsep$pe, gef.los.nsep$ci.l, gef.los.nsep$ci.u)

# Save the Results
save(t1.died.sep, t1.comp.sep , t1.los.sep,
     t1.died.nsep, t1.comp.nsep , t1.los.nsep,  
     gef.died.sep, gef.comp.sep, gef.los.sep, 
     gef.died.nsep, gef.comp.nsep, gef.los.nsep, 
     file="./routput/outcomes_est_em.RData")
     
     
     
