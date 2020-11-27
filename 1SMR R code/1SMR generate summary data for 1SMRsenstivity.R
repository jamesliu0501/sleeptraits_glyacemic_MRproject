# Overview
# This file contains the R script of generating the summary statistics for the 1SMR sensitivity analyses (i.e., collider-corrected estimates)
# Here we take sleep duration and insomnia vs Hba1c as examples  (R codes of other sleep traits are similar)
# The code of generating summary statistics of the subset SNPs of other sleep traits are similar. 
# To account for winner's curse, a process of selecting SNPs, which are significant in an independent cohort, would be added.  
# To be noticed, these estimates are in log unit which haven't be transformed into SD unit. 
# For HbA1c, you need to divide the estimates by 0.15 to obtain the SD estimates; for glucose, you need to divide the estimates by 0.17 to obtain the SD unit


library(dplyr)
library(data.table)
library(stringr)

# Continuous exposure (e.g., sleep duration)
# Summary statistics of genetic associations of with sleep duration and HbA1c in the UKB 
# (adjust for sex, age, assessment centre, chip, 40 genetic principal components)
rm(list = ls())
sleepduration <- read.table("sleepduration_hba1c_individual_data_ukb.txt", header=T)

n <- nrow(sleepduration)
size <- ncol(sleepduration[,-(1:47)]) 
G1.2 <- as.matrix(sleepduration[,-(1:47)]) 
X.2 <- as.matrix(sleepduration[,2])
Y <- as.matrix(sleepduration[,3]) 

# covariates: sex, age, assessment centre, chip, 40 genetic principal components
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



# X~G
XGdata2        = data.frame(X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   


FIT2           = summary(lm(XGdata2))                      
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= -63)                           
seBetaXG       = FIT2$coef[-1,2]                                
seBetaXG         = head(seBetaXG, n= -63)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                 


# Y~G  
YGdata        = data.frame(Y, G1.2,
                           sex, age, as.factor(centre), chip, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                           pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   

FIT3           = summary(lm(YGdata))                            
BetaYG         = FIT3$coef[-1,1]                              
BetaYG         = head(BetaYG, n= -63)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG         = head(seBetaYG, n= -63)                           


# Y~G+X 
YXGdata        = data.frame(Y, X.2, G1.2,
                            sex, age, as.factor(centre), chip, # n=4
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40) 


FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar         = head(alphahatstar, n= -63)                          
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar         = head(se.alphahatstar, n= -63)                           
betastar         = FIT4$coef[2,1] 
se.betastar  =  FIT4$coef[2,2] 


#summary statistics
# X- G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")

data<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data) <- c()

write.csv(data, "sleepduration_hba1c_adjusted_summary_data_ukb.csv", row.names = F)







# Binary exposure (e.g., insomnia)
# Summary statistics of genetic associations of with insomnia and HbA1c in the UKB 
# Genetic variant-sleep trait association was modeled using linear regression
# (adjust for sex, age, assessment centre, chip, 40 genetic principal components)

rm(list = ls())
insomnia2 <- read.table("insomnia2_hba1c_updated.txt", header=T)

n <- nrow(insomnia2)
size <- ncol(insomnia2[,-(1:47)])
G1.2 <- as.matrix(insomnia2[,-(1:47)])

# insomnia is coded as 1, 0 ;  (Never/rarely,  Sometimes) = 0,   Usually = 1
X.2 <- as.matrix(insomnia2[,2])
X.2<-ifelse(X.2<=2, 0, 1)   

Y <- as.matrix(insomnia2[,3])

# covariates: sex, age, assessment centre, chip, 40 genetic principal components
sex<-as.matrix(insomnia2[,4])
age<-as.matrix(insomnia2[,5])
centre<-as.matrix(insomnia2[,6])
chip<-as.matrix(insomnia2[,47])

pc1 <- as.matrix(insomnia2[,7])     
pc2 <- as.matrix(insomnia2[,8]) 
pc3 <- as.matrix(insomnia2[,9]) 
pc4 <- as.matrix(insomnia2[,10]) 
pc5 <- as.matrix(insomnia2[,11]) 
pc6 <- as.matrix(insomnia2[,12]) 
pc7 <- as.matrix(insomnia2[,13]) 
pc8 <- as.matrix(insomnia2[,14]) 
pc9 <- as.matrix(insomnia2[,15]) 
pc10 <- as.matrix(insomnia2[,16]) 

pc11 <- as.matrix(insomnia2[,17])     
pc12 <- as.matrix(insomnia2[,18]) 
pc13 <- as.matrix(insomnia2[,19]) 
pc14 <- as.matrix(insomnia2[,20]) 
pc15 <- as.matrix(insomnia2[,21]) 
pc16 <- as.matrix(insomnia2[,22]) 
pc17 <- as.matrix(insomnia2[,23]) 
pc18 <- as.matrix(insomnia2[,24]) 
pc19 <- as.matrix(insomnia2[,25]) 
pc20 <- as.matrix(insomnia2[,26]) 

pc21 <- as.matrix(insomnia2[,27])     
pc22 <- as.matrix(insomnia2[,28]) 
pc23 <- as.matrix(insomnia2[,29]) 
pc24 <- as.matrix(insomnia2[,30]) 
pc25 <- as.matrix(insomnia2[,31]) 
pc26 <- as.matrix(insomnia2[,32]) 
pc27 <- as.matrix(insomnia2[,33]) 
pc28 <- as.matrix(insomnia2[,34]) 
pc29 <- as.matrix(insomnia2[,35]) 
pc30 <- as.matrix(insomnia2[,36]) 

pc31 <- as.matrix(insomnia2[,37])     
pc32 <- as.matrix(insomnia2[,38]) 
pc33 <- as.matrix(insomnia2[,39]) 
pc34 <- as.matrix(insomnia2[,40]) 
pc35 <- as.matrix(insomnia2[,41]) 
pc36 <- as.matrix(insomnia2[,42]) 
pc37 <- as.matrix(insomnia2[,43]) 
pc38 <- as.matrix(insomnia2[,44]) 
pc39 <- as.matrix(insomnia2[,45]) 
pc40 <- as.matrix(insomnia2[,46]) 



# X~G
XGdata2        = data.frame(X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   

FIT2           = summary(lm(XGdata2))   # Genetic variant-sleep trait association was modeled using linear regression   
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= -63)                         
seBetaXG       = FIT2$coef[-1,2]                               
seBetaXG         = head(seBetaXG, n= -63)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                


# Y~G 
YGdata        = data.frame(Y, G1.2,
                           sex, age, as.factor(centre), chip, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                           pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   


FIT3           = summary(lm(YGdata))                          
BetaYG         = FIT3$coef[-1,1]                            
BetaYG         = head(BetaYG, n= -63)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG         = head(seBetaYG, n= -63)                          


# Y~G+X 
YXGdata        = data.frame(Y, X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40) 


FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar         = head(alphahatstar, n= -63)                           
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar         = head(se.alphahatstar, n= -63)                           
betastar         = FIT4$coef[2,1] 
se.betastar  =  FIT4$coef[2,2] 


#summary statistics
# X- G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")

data<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data) <- c()

write.csv(data, "summary_insomnia2_hba1c_ukb_ad.csv", row.names = F)








