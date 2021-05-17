# Overview
# This file contains the R script of the 1SMR sensitivity analyses (1SMRsensitivity1 (IVW), 1SMRsensitivity2 (MR-Egger), 1SMRsensitivity3 (LADreg))
# Here we take sleep duration and insomnia vs Hba1c as examples  (R codes of other sleep traits are similar)
# To be noticed, these estimates are in log unit which haven't be transformed into SD unit 
# For HbA1c, you need to divide the estimates by 0.15 log mmol/mol to obtain the SD estimates; for glucose, you need to divide the estimates by 0.17 log mmol/l to obtain the SD unit

library(foreign)
library(simex)
library(L1pack)
library(mr.raps)
library(coda)
library(xtable)

# Continuous exposure (e.g., sleep duration)
# sleep duration with HbA1c 
rm(list = ls())
sleepduration <- read.csv("summary_sleepduration_hba1c_ukb_ad.csv") 
set.seed(888)
source("ExactQ.R")  #inport "ExactQ.R" function file

# function of generating the IGx2 
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

sleepduration <- sleepduration[which(sleepduration$BetaXG>0),] 

BetaXG<-sleepduration$BetaXG 
seBetaXG<-sleepduration$seBetaXG
BetaYG<-sleepduration$BetaYG 
seBetaYG<-sleepduration$seBetaYG
alphahatstar<-sleepduration$alphahatstar 
se.alphahatstar<-sleepduration$se.alphahatstar  
betastar<-sleepduration$betastar 
se.betastar<-sleepduration$se.betastar

##############################
F = BetaXG^2/seBetaXG^2

plot(BetaXG,F)

# 1SMRsensitivity1
# Implement IVW using Modified weights 
Results = weightedIVW(BetaXG,alphahatstar,seBetaXG,se.alphahatstar,tol=0.00001)
IVW_mw  = Results$RESULTS[5,1]
Qexact  = Results$QStats[4,1]
Qp      = Results$QStats[4,2]

betaIVW = betastar[1] + IVW_mw
names(Results)


#1SMRsensitivity2
# MR-Egger   
IsqGX         = Isq(BetaXG,seBetaXG)                                                   
betahat       = summary(lm(alphahatstar~ BetaXG,weights=1/se.alphahatstar^2))$coef    
MREgger       = betastar[1] + betahat[2,1]  
# Collider correction + SiMEX 
Fit           = lm(alphahatstar~BetaXG,x=TRUE,y=TRUE,weights=1/se.alphahatstar^2)        
mod.sim2      = simex(Fit,B=500,measurement.error=seBetaXG,                              
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE")
bSIMEX        = summary(mod.sim2)$coef$jackknife   
MREggersimex  = betastar[1] + bSIMEX[2,1]
EggerSE       = bSIMEX[2,2]    
Egger_intercept = bSIMEX[1,4]


# 1SMRsensitivity3
# LAD regression
betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]   
LAD           = betastar[1] + betahat  

# Collider correction + SiMEX 
Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)         
mod.sim3      = simex(Fit,B=500,measurement.error=seBetaXG,       
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE",jackknife.estimation = FALSE)
bSIMEX        = mod.sim3$coef               
LADsimex      = betastar[1] + bSIMEX
# obtain bootstrap SE for LAD regression 
Ests = NULL
for(i in 1:100){
  L     = length(BetaXG)  
  d     = sample(L,L,replace=TRUE)  
  data3 = sleepduration[d,]   
  
  BetaXG          = data3$BetaXG
  seBetaXG        = data3$seBetaXG
  alphahatstar    = data3$alphahatstar
  se.alphahatstar = data3$se.alphahatstar
  
  # LAD regression
  betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]      
  LAD           = betastar[1] + betahat                              
  
  Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)        
  mod.sim2      = simex(Fit,B=200,measurement.error=seBetaXG,
                        SIMEXvariable="BetaXG",fitting.method="quad",
                        asymptotic="FALSE",jackknife.estimation = FALSE)
  Ests[i]       = mod.sim2$coef
  print(i)               
}

seLAD = sd(Ests) 

# sort the estimates
BetaXG          = sleepduration$BetaXG 
seBetaXG        = sleepduration$seBetaXG
alphahatstar    = sleepduration$alphahatstar
se.alphahatstar = sleepduration$se.alphahatstar

Fbar=mean(F)
Stats               = data.frame(Fbar,IsqGX,Qexact,Qp)
Estimates           = c(betaIVW,MREggersimex,LADsimex)

ColliderCorrections = Estimates-betastar[1]
SEs                 = sqrt(c(Results$RESULTS[5,2],EggerSE,seLAD)^2 + (se.betastar[1])^2)
LCIs                = Estimates - 1.96 * SEs
UCIs                = Estimates + 1.96 * SEs
pval                = 2*(1-pnorm(abs(Estimates/SEs)))
Final = data.frame(Estimates,SEs,pval,
                   row.names=c("IVW","MR-Egger","LADreg"))
fsNigx<-c(Fbar,IsqGX,NA)

sleepduration_hba1c = data.frame(Estimates, LCIs, UCIs, pval,fsNigx,
                                row.names=c("IVW","MR-Egger","LADreg"))






# Binary exposure (e.g., insomnia)
# Insomnia with HbA1c 
rm(list = ls())
insomnia2 <- read.csv("summary_insomnia2_hba1c_ukb_ad.csv")

set.seed(888)
source("ExactQ.R") #inport "ExactQ.R" function file
# function of generating the IGx2 
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

insomnia2 <- insomnia2[which(insomnia2$BetaXG>0),] 
BetaXG<-insomnia2$BetaXG
seBetaXG<-insomnia2$seBetaXG
BetaYG<-insomnia2$BetaYG
seBetaYG<-insomnia2$seBetaYG
alphahatstar<-insomnia2$alphahatstar
se.alphahatstar<-insomnia2$se.alphahatstar  
betastar<-insomnia2$betastar
se.betastar<-insomnia2$se.betastar

F = BetaXG^2/seBetaXG^2
plot(BetaXG,F)

# 1SMRsensitivity1
# Implement IVW using Modified weights  
Results = weightedIVW(BetaXG,alphahatstar,seBetaXG,se.alphahatstar,tol=0.00001)
IVW_mw  = Results$RESULTS[5,1]
Qexact  = Results$QStats[4,1]
Qp      = Results$QStats[4,2]
betaIVW = betastar[1] + IVW_mw
names(Results)


# 1SMRsensitivity2
# MR-Egger 
IsqGX         = Isq(BetaXG,seBetaXG)                                                  
betahat       = summary(lm(alphahatstar~ BetaXG,weights=1/se.alphahatstar^2))$coef   
MREgger       = betastar[1] + betahat[2,1]  
#Collider correction + SiMEX 
Fit           = lm(alphahatstar~BetaXG,x=TRUE,y=TRUE,weights=1/se.alphahatstar^2)        
mod.sim2      = simex(Fit,B=500,measurement.error=seBetaXG,                             
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE")
bSIMEX        = summary(mod.sim2)$coef$jackknife   
MREggersimex  = betastar[1] + bSIMEX[2,1]
EggerSE       = bSIMEX[2,2]    


# 1SMRsensitivity3
# LAD regression
betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]   
LAD           = betastar[1] + betahat                            

Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)         
mod.sim3      = simex(Fit,B=500,measurement.error=seBetaXG,      
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE",jackknife.estimation = FALSE)
bSIMEX        = mod.sim3$coef               
LADsimex      = betastar[1] + bSIMEX

# obtain bootstrap SE for LAD regression 
Ests = NULL
for(i in 1:100){
  L     = length(BetaXG)  
  d     = sample(L,L,replace=TRUE)  
  data3 = insomnia2[d,]   
  
  BetaXG          = data3$BetaXG
  seBetaXG        = data3$seBetaXG
  alphahatstar    = data3$alphahatstar
  se.alphahatstar = data3$se.alphahatstar
  
  # LAD regression
  betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]      
  LAD           = betastar[1] + betahat                               
  
  Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)        
  mod.sim2      = simex(Fit,B=200,measurement.error=seBetaXG,
                        SIMEXvariable="BetaXG",fitting.method="quad",
                        asymptotic="FALSE",jackknife.estimation = FALSE)
  Ests[i]       = mod.sim2$coef
  print(i)               
}

seLAD = sd(Ests) 


# sort the estimates
BetaXG          = insomnia2$BetaXG 
seBetaXG        = insomnia2$seBetaXG
alphahatstar    = insomnia2$alphahatstar
se.alphahatstar = insomnia2$se.alphahatstar


Fbar=mean(F)
Stats               = data.frame(Fbar,IsqGX,Qexact,Qp)
Estimates           = c(betaIVW,MREggersimex,LADsimex)
ColliderCorrections = Estimates-betastar[1]
SEs                 = sqrt(c(Results$RESULTS[5,2],EggerSE,seLAD)^2 + (se.betastar[1])^2)
LCIs                = Estimates - 1.96 * SEs
UCIs                = Estimates + 1.96 * SEs
pval                = 2*(1-pnorm(abs(Estimates/SEs)))

Final = data.frame(Estimates,SEs,pval,
                   row.names=c("IVW","MR-Egger","LADreg"))

fsNigx<-c(Fbar,IsqGX,NA)

insomnia2_hba1c_all = data.frame(Estimates, LCIs, UCIs, pval,fsNigx,
                                 row.names=c("IVW","MR-Egger","LADreg"))

