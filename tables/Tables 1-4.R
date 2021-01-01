

rm(list=ls())

source("table_printing.R")

load("outcomes_est_em.RData")

# Table 1
load("iv_balance_match.RData")

out.tb.post                               
xtable(out.tb.post)


# Table 2
# IV Estimates
print_out_table(gef.died.sep*100, gef.comp.sep*100, gef.los.sep, 
                gef.died.nsep*100, gef.comp.nsep*100, gef.los.nsep)
   
# Table 3
rm(list=ls())

load("sens_results.RData")

cat("Complication  &", round(sen.comp.sep$Changepoint, 2), "&",  round(sen.comp.nsep$Changepoint, 2), " \\\\
     Length of Stay  &", round(sen.los.sep$Changepoint, 2), "&",  round(sen.los.nsep$Changepoint, 2),  "\\\\ \n")
    
# Table 4                  
rm(list=ls())

load("sens_results.RData")
load("sens_results_effmod.RData")


# LOS
# Sepsis Group
cat("Conventional & [", round(sen.los.sep$SensitivityInterval[,2], 2), ",", round(sen.los.sep$SensitivityInterval[,3], 2), "] & ",  
        round(sen.los.sep$Changepoint, 2), "\\\\ \n")
cat("Regression-based & [", round(sen.los.sep.lin.het$SensitivityInterval[,2], 2), ",", round(sen.los.sep.lin.het$SensitivityInterval[,3], 2), "] & ",  
        round(sen.los.sep.lin.het$Changepoint, 2), "\\\\ \n")
cat("Pair of Pairs SE & [", round(sen.los.sep.het$SensitivityInterval[,2], 2), ",", round(sen.los.sep.het$SensitivityInterval[,3], 2), "] & ",  
        round(sen.los.sep.het$Changepoint, 2), "\\\\ \n")

# Non-Sepsis Group
cat("Conventional  &  [", round(sen.los.nsep$SensitivityInterval[,2], 2), ",", round(sen.los.nsep$SensitivityInterval[,3], 2), "] & ",  
        round(sen.los.nsep$Changepoint, 2), "\\\\ \n")
cat("Regression-based  & [", round(sen.los.nsep.lin.het$SensitivityInterval[,2], 2), ",", round(sen.los.nsep.lin.het$SensitivityInterval[,3], 2), "] & ",  
        round(sen.los.nsep.lin.het$Changepoint, 2), "\\\\ \n")
cat("air of Pairs SE  & [", round(sen.los.nsep.het$SensitivityInterval[,2], 2), ",", round(sen.los.nsep.het$SensitivityInterval[,3], 2), "] & ",  
        round(sen.los.nsep.het$Changepoint, 2), "\\\\ \n")


# p-values reported in the Sec 6.2
sen.los.sep.het.pval
sen.los.nsep.het.pval

