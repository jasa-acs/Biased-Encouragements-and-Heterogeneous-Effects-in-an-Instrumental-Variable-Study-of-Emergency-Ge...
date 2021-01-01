# Biased Encouragements and Heterogeneous Effects in an Instrumental Variable Study of Emergency General Surgical Outcomes

# Author Contributions Checklist Form

## Data

### Abstract 

The data link the American Medical Association (AMA) Physician Masterfile with all-payer hospital discharge claims from New York, Florida and Pennsylvania in 2012-2013. The study population consists of all patients admitted for inpatient care emergently, urgently, or through the emergency department. The data also include patient sociodemographic and clinical characteristics including indicators for frailty, severe sepsis  or septic shock, and 31 comorbidities based on Elixhauser indices.

###Availability

The data used in this paper are not publicly available because they contain sensitive information about patients’ medical history. They can, however, be accessed from each state through a formal application.  We include a simulated data set to allow users run the functions without access to the original data used in our analysis.

## Code

### Description 

The R code to run all the analyses in the manuscript along with supporting functions for are available in GitHub as one repository.  The repository contains a large number of R functions that performs specific analyses. R code for all simulations in the paper are also included in the repository.

### Optional Information 

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

## Instructions for use

### Reproducibility 

A README file is included.  It describes the scripts used for the analysis in the manuscript.  It also describes a more succinct analysis based on included simulated data.  This analysis allows users to estimate IV effects and perform the key types of sensitivity analysis described in the paper.  It also describes two scripts that reproduce all the simulations reported in the paper.
