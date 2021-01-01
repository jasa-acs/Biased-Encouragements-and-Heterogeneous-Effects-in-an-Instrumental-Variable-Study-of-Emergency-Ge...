
library(foreign)
library(xtable)

rm(list=ls())

# Load Data
data <- read.dta("./data/iv_match_s1.dta")
data <- data[order(data$pref_bin, decreasing = TRUE), ]

# Define Quantities for Balance Checking
data$id <- seq(1:nrow(data))
data$treat <- data$pref_bin
n.pat <- table(data$pref_bin)

# Load R Functions
load("./rfunctions/iv_match_final_em.RData")
source("./rfunctions/meantabRC.R")

bal_vars <- c("id", "treat","age","comorb", "agenow", "yearsexp",
              "female", "hisp", "white", "afam", "other", 
             "angus", "disability_bin", 
             "ynel1", "ynel2","ynel3","ynel4", "ynel5","ynel6", "ynel7",
             "ynel8","ynel9", "ynel10","ynel11","ynel12","ynel13","ynel14",
             "ynel15","ynel16","ynel17", "ynel18","ynel19","ynel20","ynel21",
             "ynel22","ynel23","ynel24","ynel25","ynel26", "ynel27","ynel28",
             "ynel29","ynel30","ynel31")

# Calc Balance on Unmatched Data  
X_mat = data[bal_vars]
ids <- data$id
bal.tab.pre <- meantabRC(X_mat, ids, exact=NULL, digits=2)
bal.tab.pre

# Calc Balance on Matched Data
B_mat = match.data[bal_vars]
ids <- match.data$id
bal.tab.post <- meantabRC(B_mat, ids, exact=NULL, digits=2)
bal.tab.post

# Collect Results into One Table
out.tb.pre <- bal.tab.pre[,4:6]
out.tb.post <- bal.tab.post[,4:6]
out.tb <- cbind(out.tb.pre, out.tb.post)
out.tb

# Variable Names
var.names <- c("Age", "No. Comorbidities", "Surgeon Age", "Surgeon Experience (Years)",
               "Female", "Hispanic", "White", "African-American", 
               "Other Racial Cat.", "Sepsis", "Disability", 
               "Congestive Heart Failure",  "Cardiac Arrhythmias","Valvular Disease",
               "Pulmonary Circulation Disorders", "Peripheral Vascular Disorders",
               "Hypertension, Uncomplicated","Paralysis", "Other Neurological Disorders", 
               "Chronic Pulmonary Disease","Diabetes, Uncomplicated",
               "Diabetes, Complicated","Hypothyroidism", "Renal Failure","Liver Disease",
               "Peptic Ulcer Disease Excluding Bleeding","AIDS/HIV", "Lymphoma","Metastatic Cancer",
               "Solid Tumor Without Metastasis","Rheumatoid Arthritis/Collagen Vascular",
               "Coagulopathy","Obesity","Weight Loss","Fluid and Electrolyte Disorders",
               "Blood Loss Anemia","Deficiency Anemia","Alcohol Abuse","Drug Abuse",
               "Psychoses","Depression","Hypertension, Complicated")
               
rownames(out.tb) <- var.names
colnames(out.tb) <- c("High Operative Pref", "Low Operative Pref", "Standardized Difference",
                      "High Operative Pref", "Low Operative Pref", "Standardized Difference")                            
out.tb.post <- out.tb                

# Sample Sizes Before and After
n.pat
n.pat.post <- table(match.data$treat)
n.pat.post

save(out.tb.post, n.pat.post, n.pat, file ="./routput/iv_balance_match.RData")

