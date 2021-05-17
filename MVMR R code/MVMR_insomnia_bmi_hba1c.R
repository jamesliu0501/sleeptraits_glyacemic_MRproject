# Overview
# This file contains the R script of the sensitivity multivariable Mendelian randomization (MVMR) of insomnia and BMI on HbA1c
# To be noticed, these estimates are in log unit which haven't be transformed into SD unit. 
# For HbA1c, you need to divide the estimates by 0.15 log mmol/mol to obtain the estimates in SD unit

library(dplyr)
library(data.table)
library(stringi)
library(ff)
library(AER)
library(ivpack)

rm(list = ls())

data<- fread(file="data.csv", fill=T, header = T)

# continuous covariates (baseline age + PC1 - PC40)
rc_numeric<-cbind(data$age_recruitment,
                  data$PC01,data$PC02,data$PC03,data$PC04,data$PC05,data$PC06,data$PC07,data$PC08,data$PC09,data$PC010,
                  data$PC011,data$PC012,data$PC013,data$PC014,data$PC015,data$PC016,data$PC017,data$PC018,data$PC019,data$PC020,
                  data$PC021,data$PC022,data$PC023,data$PC024,data$PC025,data$PC026,data$PC027,data$PC028,data$PC029,data$PC030,
                  data$PC031,data$PC032,data$PC033,data$PC034,data$PC035,data$PC036,data$PC037,data$PC038,data$PC039,data$PC040)



# 1st stage to predict insomnia 
# (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
data$preinsomnia <- predict(glm(insomnia_binary ~ insomnia_unweighted2 + bmi_unweighted + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,
                                na.action = na.exclude, data=data))

# 1st stage to predict bmi 
# (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
data$prebmi <- predict(glm(bmi ~ insomnia_unweighted2 + bmi_unweighted + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,
                           na.action = na.exclude, data=data))

# 2nd stage
# (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
twoStage <- glm(log_hba1c ~ preinsomnia + prebmi  + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, 
                data=data)

insomnia_bmi_hba1c<-coef(summary(twoStage))
n_insomnia_bmi_hba1c_mvmr<-summary(twoStage)$'df'[1] + summary(twoStage)$'df'[2]

# Insomnia
insomnia_hba1c_mvmr_b<-insomnia_bmi_hba1c[2,1]
insomnia_hba1c_mvmr_lci<-insomnia_bmi_hba1c[2,1] - 1.96 * insomnia_bmi_hba1c[2,2]
insomnia_hba1c_mvmr_uci<-insomnia_bmi_hba1c[2,1] + 1.96 * insomnia_bmi_hba1c[2,2]
insomnia_hba1c_mvmr_p<-insomnia_bmi_hba1c[2,4]

# BMI
bmi_hba1c_mvmr_b<-insomnia_bmi_hba1c[3,1]
bmi_hba1c_mvmr_lci<-insomnia_bmi_hba1c[3,1] - 1.96 * insomnia_bmi_hba1c[3,2]
bmi_hba1c_mvmr_uci<-insomnia_bmi_hba1c[3,1] + 1.96 * insomnia_bmi_hba1c[3,2]
bmi_hba1c_mvmr_p<-insomnia_bmi_hba1c[3,4]





# Univariable Mendelian randomization (UVMR) of BMI with HbA1c
# BMI only ( adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
mr<-ivreg(log_hba1c ~ bmi + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric 
          | bmi_unweighted + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,  
          data = data)

n_bmi_hba1c<-summary(mr)$'df'[1] + summary(mr)$'df'[2]
bmi_hba1c<-coef(summary(mr))
bmi_hba1c_2sls_b<-bmi_hba1c[2,1]
bmi_hba1c_2sls_lci<-bmi_hba1c[2,1] - 1.96 * bmi_hba1c[2,2]
bmi_hba1c_2sls_uci<-bmi_hba1c[2,1] + 1.96 * bmi_hba1c[2,2]
bmi_hba1c_2sls_p<-bmi_hba1c[2,4]






# Sanderson–Windmeijer conditional F statistics of the unweighted allele score of insomnia and BMI https://doi.org/10.1093/ije/dyy262
rm(list = ls())
data_updated<-read.csv("MVMR conditional f statistics.csv") # (A complete dataset, excluding missing data)

rc_numeric<-cbind(data_updated$age_recruitment,
                  data_updated$PC01,data_updated$PC02,data_updated$PC03,data_updated$PC04,data_updated$PC05,data_updated$PC06,data_updated$PC07,data_updated$PC08,data_updated$PC09,data_updated$PC010,
                  data_updated$PC011,data_updated$PC012,data_updated$PC013,data_updated$PC014,data_updated$PC015,data_updated$PC016,data_updated$PC017,data_updated$PC018,data_updated$PC019,data_updated$PC020,
                  data_updated$PC021,data_updated$PC022,data_updated$PC023,data_updated$PC024,data_updated$PC025,data_updated$PC026,data_updated$PC027,data_updated$PC028,data_updated$PC029,data_updated$PC030,
                  data_updated$PC031,data_updated$PC032,data_updated$PC033,data_updated$PC034,data_updated$PC035,data_updated$PC036,data_updated$PC037,data_updated$PC038,data_updated$PC039,data_updated$PC040)


# conditional F statistics for insomnia
reg1 <- ivreg(insomnia_binary ~ bmi + rc_numeric|bmi_unweighted + insomnia_unweighted2 + rc_numeric, data = data_updated)
residuals <- residuals(reg1)

res1 <-(lm(residuals ~ rc_numeric))
res2 <- (lm(residuals ~ bmi_unweighted + insomnia_unweighted2 + rc_numeric, data = data_updated))

F_insomnia <- anova(res1, res2)$F[2]*(anova(res1,res2)$Df[2]/(anova(res1,res2)$Df[2]-1))




# conditional F statistics for BMI
reg1 <- ivreg(bmi ~ insomnia_binary + rc_numeric|bmi_unweighted + insomnia_unweighted2 + rc_numeric, data = data_updated)
residuals <- residuals(reg1)

res1 <-(lm(residuals ~ rc_numeric))
res2 <- (lm(residuals ~ bmi_unweighted + insomnia_unweighted2 + rc_numeric, data = data_updated))

F_bmi <- anova(res1, res2)$F[2]*(anova(res1,res2)$Df[2]/(anova(res1,res2)$Df[2]-1))





