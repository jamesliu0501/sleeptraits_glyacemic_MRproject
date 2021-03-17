# Overview
# This file contains the R script of the observational multivariable regression (MV) of Model 1 and Model 2
# Here we take sleep duration and insomnia vs HbA1c as examples  (R codes of other sleep traits are similar)
# For glucose, fasting time and dilution factor would be additionally adjusted for 
# To be noticed, these estimates are in log unit which haven't be transformed into SD unit. 
# For HbA1c, you need to divide the estimates by 0.15 to obtain the SD estimates; for glucose, you need to divide the estimates by 0.17 to obtain the SD unit

library(dplyr)
library(data.table)
library(stringi)
library(ff)
library(olsrr)

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
# sleep duration vs HbA1c
# Model 1: adjusted for sex, age, assessment centre, 40 genetic principal components, smoking, alcohol intake, Townsend deprivation, education, and physical activity

mvr_a_log<-glm(log_hba1c~ sleepduration + 
                 sex + as.factor(assessment_centre) +
                 rc_numeric +  # age and PC1-40 are included here
                 townsend + 
                 as.factor(smoking) + as.factor(aif) + as.factor(degree) + as.factor(vig),
               data = data)


sleepduration_a_log_hba1c<-coef(summary(mvr_a_log))

sleepduration_a_log_hba1c_b<-sleepduration_a_log_hba1c[2,1]
sleepduration_a_log_hba1c_lci<-sleepduration_a_log_hba1c[2,1] - 1.96 * sleepduration_a_log_hba1c[2,2]
sleepduration_a_log_hba1c_uci<-sleepduration_a_log_hba1c[2,1] + 1.96 * sleepduration_a_log_hba1c[2,2]
sleepduration_a_log_hba1c_p<-sleepduration_a_log_hba1c[2,4]

n_sleepduration_a_log_hba1c<-summary(mvr_a_log)$'df.null'+1


# Model 2: additionally adjusted for BMI.
mvr_a2_log<-glm(log_hba1c~ sleepduration + 
                  sex + as.factor(assessment_centre) +
                  rc_numeric + # age and PC1-40 are included here
                  townsend + bmi +
                  as.factor(smoking) + as.factor(aif) + as.factor(degree) + as.factor(vig) ,
                data = data)


sleepduration_a2_log_hba1c<-coef(summary(mvr_a2_log))

sleepduration_a2_log_hba1c_b<-sleepduration_a2_log_hba1c[2,1]
sleepduration_a2_log_hba1c_lci<-sleepduration_a2_log_hba1c[2,1] - 1.96 * sleepduration_a2_log_hba1c[2,2]
sleepduration_a2_log_hba1c_uci<-sleepduration_a2_log_hba1c[2,1] + 1.96 * sleepduration_a2_log_hba1c[2,2]
sleepduration_a2_log_hba1c_p<-sleepduration_a2_log_hba1c[2,4]

n_sleepduration_a2_log_hba1c<-summary(mvr_a2_log)$'df.null'+1


# Binary exposure (e.g., insomnia)
# insomnia vs HbA1c
# Model 1: adjusted for sex, age, assessment centre, 40 genetic principal components, smoking, alcohol intake, Townsend deprivation, education, and physical activity

mvr_a_log<-glm(log_hba1c~ insomnia_binary + 
                 sex + as.factor(assessment_centre) +
                 rc_numeric + 
                 townsend + 
                 as.factor(smoking) + as.factor(aif) + as.factor(degree) + as.factor(vig),
               data = data)

insomnia_a_log_hba1c<-coef(summary(mvr_a_log))

insomnia_a_log_hba1c_b<-insomnia_a_log_hba1c[2,1]
insomnia_a_log_hba1c_lci<-insomnia_a_log_hba1c[2,1] - 1.96 * insomnia_a_log_hba1c[2,2]
insomnia_a_log_hba1c_uci<-insomnia_a_log_hba1c[2,1] + 1.96 * insomnia_a_log_hba1c[2,2]
insomnia_a_log_hba1c_p<-insomnia_a_log_hba1c[2,4]

n_insomnia_a_log_hba1c<-summary(mvr_a_log)$'df.null'+1


# Model 2: additionally adjusted for BMI.
mvr_a2_log<-glm(log_hba1c~ insomnia_binary + 
                  sex + as.factor(assessment_centre) +
                  rc_numeric + 
                  townsend + bmi +
                  as.factor(smoking) + as.factor(aif) + as.factor(degree) + as.factor(vig),
                data = data)

insomnia_a2_log_hba1c<-coef(summary(mvr_a2_log))

insomnia_a2_log_hba1c_b<-insomnia_a2_log_hba1c[2,1]
insomnia_a2_log_hba1c_lci<-insomnia_a2_log_hba1c[2,1] - 1.96 * insomnia_a2_log_hba1c[2,2]
insomnia_a2_log_hba1c_uci<-insomnia_a2_log_hba1c[2,1] + 1.96 * insomnia_a2_log_hba1c[2,2]
insomnia_a2_log_hba1c_p<-insomnia_a2_log_hba1c[2,4]

n_insomnia_a2_log_hba1c<-summary(mvr_a2_log)$'df.null'+1






# Exclusion 1 
## Exclude participants with diabetic status defined by the Eastwood algorithm
rm(list = ls())

data<- fread(file="data.csv",
             fill=T, header = T)

data<-data.table(data)
data<-subset(data, adj_probposs_t1dm == 0 & adj_probposs_t2dm == 0)


# Exclusion 2 
## Additionally exclude participants with HbA1c >= 48 mmol/mol
#exclusion (exclude participants with self-reported diabetes or insulin use at the baseline)
# same R code could be applied as above 
rm(list = ls())
data<- fread(file="data.csv",
             fill=T, header = T)

data<-data.table(data)
data<-subset(data, adj_probposs_t1dm == 0 & adj_probposs_t2dm == 0)
data<-data[data$hba1c1<48,] 


