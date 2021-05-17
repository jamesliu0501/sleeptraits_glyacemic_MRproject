# Overview
# This file contains the R script of the Sargan test of the 1SMRmain (2SLS) analyses to explore between SNP heterogeneity 
# Here we take sleep duration vs HbA1c as an example  (R codes of other sleep traits are similar)

library(dplyr)
library(data.table)
library(stringi)
library(ff)
library(AER)

rm(list = ls())

sleepduration <- read.table("sleepduration_hba1c_individual_data_ukb.txt", header=T)

#####
G1.2 <- as.matrix(sleepduration[,-(1:47)])
X.2 <- as.matrix(sleepduration[,2])
Y <- as.matrix(sleepduration[,3])

# Covariates: sex, age, assessment centre, chip, 40 genetic principal components
sex<-as.matrix(sleepduration[,4])
age<-as.matrix(sleepduration[,5])
centre<-as.matrix(sleepduration[,6])
chip<-as.matrix(sleepduration[,47])

pc1 <- as.matrix(sleepduration[,7])     
pc2 <- as.matrix(sleepduration[,8]) 
pc3 <- as.matrix(sleepduration[,9]) 
pc4 <- as.matrix(sleepduration[,10]) 
pc5 <- as.matrix(sleepduration[,11]) 
pc6 <- as.matrix(sleepduration[,12]) 
pc7 <- as.matrix(sleepduration[,13]) 
pc8 <- as.matrix(sleepduration[,14]) 
pc9 <- as.matrix(sleepduration[,15]) 
pc10 <- as.matrix(sleepduration[,16]) 

pc11 <- as.matrix(sleepduration[,17])     
pc12 <- as.matrix(sleepduration[,18]) 
pc13 <- as.matrix(sleepduration[,19]) 
pc14 <- as.matrix(sleepduration[,20]) 
pc15 <- as.matrix(sleepduration[,21]) 
pc16 <- as.matrix(sleepduration[,22]) 
pc17 <- as.matrix(sleepduration[,23]) 
pc18 <- as.matrix(sleepduration[,24]) 
pc19 <- as.matrix(sleepduration[,25]) 
pc20 <- as.matrix(sleepduration[,26]) 

pc21 <- as.matrix(sleepduration[,27])     
pc22 <- as.matrix(sleepduration[,28]) 
pc23 <- as.matrix(sleepduration[,29]) 
pc24 <- as.matrix(sleepduration[,30]) 
pc25 <- as.matrix(sleepduration[,31]) 
pc26 <- as.matrix(sleepduration[,32]) 
pc27 <- as.matrix(sleepduration[,33]) 
pc28 <- as.matrix(sleepduration[,34]) 
pc29 <- as.matrix(sleepduration[,35]) 
pc30 <- as.matrix(sleepduration[,36]) 

pc31 <- as.matrix(sleepduration[,37])     
pc32 <- as.matrix(sleepduration[,38]) 
pc33 <- as.matrix(sleepduration[,39]) 
pc34 <- as.matrix(sleepduration[,40]) 
pc35 <- as.matrix(sleepduration[,41]) 
pc36 <- as.matrix(sleepduration[,42]) 
pc37 <- as.matrix(sleepduration[,43]) 
pc38 <- as.matrix(sleepduration[,44]) 
pc39 <- as.matrix(sleepduration[,45]) 
pc40 <- as.matrix(sleepduration[,46]) 

## Adjust for sex, age, assessment centre, chip, 40 genetic principal components
mr<-ivreg(Y ~ X.2 + sex + age + as.factor(centre) + chip + 
            pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 +
            pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 +
            pc21 + pc22 + pc23 + pc24 + pc25 + pc26 + pc27 + pc28 + pc29 + pc30 +
            pc31 + pc32 + pc33 + pc34 + pc35 + pc36 + pc37 + pc38 + pc39 + pc40   
          | G1.2 + sex + age + as.factor(centre) + chip + 
            pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 +
            pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 +
            pc21 + pc22 + pc23 + pc24 + pc25 + pc26 + pc27 + pc28 + pc29 + pc30 +
            pc31 + pc32 + pc33 + pc34 + pc35 + pc36 + pc37 + pc38 + pc39 + pc40 )

summary<-summary(mr, diagnostics=TRUE)

coef<-summary$coefficients

sleepduration_hba1c_b<-coef[2,1]
sleepduration_hba1c_lci<-coef[2,1] - 1.96 * coef[2,2]
sleepduration_hba1c_uci<-coef[2,1] + 1.96 * coef[2,2]
sleepduration_hba1c_p<-coef[2,4]

sleepduration_hba1c_n_snp<-summary$diagnostics[3,1] + 1
sleepduration_hba1c_sargan_p<-summary$diagnostics[3,4]



