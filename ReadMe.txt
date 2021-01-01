

This replication file includes all the computer code used to generate the results in the paper. We are, however, unable to share either the original data or the summaries generated from the analysis. However, we include a script that uses simulated data to conduct all the main forms of analysis in the paper.

Overview and Description of Directories and Scripts

./code
  /Application - R scripts for the analysis conducted in the paper. We provide notes on those scripts below. The order of the scrips follows the order 
required for the analysis.
	

-IV_Sens_Match_ER.R -- In this script, we conduct the near-far match and save the matched data. We exactly matched on hospitals and sepsis status.
-Balance_Statistics.r -- In this script, we measured balance before and after the near-far match.
-IV_Sens_Outcome.R -- In this script, we estimate the ITT and IV estimates after matching. We save point estimates and confidence intervals.
-IV_Sens.R -- In this script, we conduct the sensitivity analyses that do not account for possible effect modification.
-IV_Sens_Eff_Mod.R -- In this script, we conduct the sensitivity analyses that account for possible effect modification.

/Analysis - Example Data 

-IV_Sens_Example.R -- In this script, we conduct a self-contained analysis based on simulated data. The simulated data has a matched pair structure.
In the script, we estimate the IV effect using the effect ratio, conduct a sensitivity analysis that does not account for possible effect modification,
and a second sensitivity analysis that accounts for effect modification. It uses the same functions used in the main analysis 

/Simulation - R scripts for both simulations in the paper. These are labeled by which table is produced in the paper.

./data - contains simulated data only.  Does not contain actual data used in the analysis.

Description of Simulated Data: ExampleData.csv

Sample Size - 5000
Z - Binary IV 
D - Binary Treatment Variable
Y - Outcome
X.1 - X.10 - Baseline Covariates Drawn from Uniform Distributions.
Pair - Matched Pair IDs

./rfunctions - all R functions used in the analysis
./routput - Empty directory, we are unable to share analysis output
./tables - R script which generates the tables in the paper.

Computing Information:
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

attached packages:
rcbalance_1.8.5 plyr_1.8.6      
MASS_7.3-51.6   
optmatch_0.9-13 
survival_3.1-12 
xtable_1.8-4    
foreign_0.8-80  
dplyr_1.0.0    






