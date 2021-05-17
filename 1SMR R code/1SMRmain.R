# Overview
# This file contains the R script of the main analysis of one-sample Mendelian randomization (1SMRmain) (2SLS)
# Here we take sleep duration and insomnia vs HbA1c as examples (R codes of other sleep traits are similar)
# For glucose, fasting time and dilution factor would be additionally adjusted for 
# To be noticed, these estimates are in log unit which haven't be transformed into the SD unit. 
# For HbA1c, you need to divide the estimates by 0.15 log mmol/mol to obtain the SD estimates; for glucose, you need to divide the estimates by 0.17 log mmol/lto obtain the SD unit

library(dplyr)
library(data.table)
library(stringi)
library(ff)
library(AER)
library(ivpack)

# Read the UKB data
data<- fread(file="data.csv",
             fill=T, header = T)

#continuous covariates (baseline age + PC1 - PC40)
rc_numeric<-cbind(data$age_recruitment,
                  data$PC01,data$PC02,data$PC03,data$PC04,data$PC05,data$PC06,data$PC07,data$PC08,data$PC09,data$PC010,
                  data$PC011,data$PC012,data$PC013,data$PC014,data$PC015,data$PC016,data$PC017,data$PC018,data$PC019,data$PC020,
                  data$PC021,data$PC022,data$PC023,data$PC024,data$PC025,data$PC026,data$PC027,data$PC028,data$PC029,data$PC030,
                  data$PC031,data$PC032,data$PC033,data$PC034,data$PC035,data$PC036,data$PC037,data$PC038,data$PC039,data$PC040)


# Continuous exposure (e.g., sleep duration)
# Sleep duration vs HbA1c (log unit)
# 1SMRmain using 2sls (1SMRmain adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
mr<-ivreg(log_hba1c ~ sleepduration  + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric 
          | sleepduration_unweighted + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,  
          data = data)

sleepduration_hba1c<-coef(summary(mr))

sleepduration_hba1c_2sls_b<-sleepduration_hba1c[2,1]
sleepduration_hba1c_2sls_lci<-sleepduration_hba1c[2,1] - 1.96 * sleepduration_hba1c[2,2]
sleepduration_hba1c_2sls_uci<-sleepduration_hba1c[2,1] + 1.96 * sleepduration_hba1c[2,2]
sleepduration_hba1c_2sls_p<-sleepduration_hba1c[2,4]

n_sleepduration_hba1c<-summary(mr)$'df'[1] + summary(mr)$'df'[2]




# Binary exposure (e.g., insomnia)
# Insomnia vs HbA1c (log unit)
# linear regression of IV - insomnia: representing the average difference in HbA1c if everyone suffer from "usually" insomnia compared to if everyone suffer from "sometimes" or "rarely/never" insomnia 
# 1SMRmain using 2sls (1SMRmain adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
mr<-ivreg(log_hba1c ~ insomnia_binary + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric 
          | insomnia_unweighted2 + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,  
          data = data)

insomnia2_hba1c<-coef(summary(mr))

insomnia2_hba1c_2sls_b<-insomnia2_hba1c[2,1]
insomnia2_hba1c_2sls_lci<-insomnia2_hba1c[2,1] - 1.96 * insomnia2_hba1c[2,2]
insomnia2_hba1c_2sls_uci<-insomnia2_hba1c[2,1] + 1.96 * insomnia2_hba1c[2,2]
insomnia2_hba1c_2sls_p<-insomnia2_hba1c[2,4]

n_insomnia2_hba1c<-summary(mr)$'df'[1] + summary(mr)$'df'[2]







# Exclusion 1 
## Exclude participants with with diabetic status defined by the Eastwood algorithm
rm(list = ls())

data<- fread(file="data.csv",
             fill=T, header = T)

data<-data.table(data)
data<-subset(data, adj_probposs_t1dm == 0 & adj_probposs_t2dm == 0)

# Exclusion 2 
## Additionally exclude participants with HbA1c >= 48 mmol/mol
# same R code could be applied as above 
rm(list = ls())
data<- fread(file="data.csv",
             fill=T, header = T)

data<-data.table(data)
data<-subset(data, adj_probposs_t1dm == 0 & adj_probposs_t2dm == 0)
data<-data[data$hba1c1<48,] 


# same R code could be applied as above 





