# Overview
# This file contains the R script of the sensitivity one-sample Mendelian randomization of HbA1c on insomnia

library(dplyr)
library(data.table)
library(stringi)
library(ff)
library(AER)

rm(list = ls())
data<- fread(file="data.csv",
             fill=T, header = T)

# continuous covariates (baseline age + PC1 - PC40)
rc_numeric<-cbind(data$age_recruitment,
                  data$PC01,data$PC02,data$PC03,data$PC04,data$PC05,data$PC06,data$PC07,data$PC08,data$PC09,data$PC010,
                  data$PC011,data$PC012,data$PC013,data$PC014,data$PC015,data$PC016,data$PC017,data$PC018,data$PC019,data$PC020,
                  data$PC021,data$PC022,data$PC023,data$PC024,data$PC025,data$PC026,data$PC027,data$PC028,data$PC029,data$PC030,
                  data$PC031,data$PC032,data$PC033,data$PC034,data$PC035,data$PC036,data$PC037,data$PC038,data$PC039,data$PC040)

data$insomnia_binary_updated<-data$insomnia_binary - 1  # to code insomnia_binary as 0, 1; (Never/rarely,  Sometimes) = 0,   Usually = 1

# 1st stage (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
data$prehba1c <- predict(glm(log_hba1c ~ hba1c_unweighted + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,
                             na.action = na.exclude, data=data))

# 2nd stage (logistic regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
twoStage <- glm(insomnia_binary_updated ~ prehba1c + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, 
                family=binomial(link='logit'), data=data)

hba1c_insomnia<-coef(summary(twoStage))


hba1c_insomnia_b<-exp(hba1c_insomnia[2,1])
hba1c_insomnia_lci<-exp(hba1c_insomnia[2,1] - 1.96 * hba1c_insomnia[2,2])
hba1c_insomnia_uci<-exp(hba1c_insomnia[2,1] + 1.96 * hba1c_insomnia[2,2])
hba1c_insomnia_p<-hba1c_insomnia[2,4]

n_hba1c_insomnia<-summary(twoStage)$'df'[1] + summary(twoStage)$'df'[2]

