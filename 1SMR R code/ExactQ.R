weightedIVW = function(BetaXG,BetaYG,seBetaXG,seBetaYG,tol=0.00001){

F             = BetaXG^2/seBetaXG^2
mF            = mean(F)
DF            = length(BetaYG)-1
BIV           = BetaYG/BetaXG


IVWFunc = function(W=1){

BIVw      = BIV*sqrt(W)
sW        = sqrt(W)
IVWfit    = summary(lm(BIVw ~ -1+sW))
phi_IVW   = IVWfit$sigma^2

# Inference with correct standard errors

Bhat    = IVWfit$coef[1,1]
SE      = IVWfit$coef[1,2]/min(1,phi_IVW)
IVW_p   = 2*(1-pt(abs(Bhat/SE),DF))
IVW_CI  = Bhat + c(-1,1)*qt(df=DF, 0.975)*SE

# Q statistics

QIVW     = DF*phi_IVW
Qp       = 1-pchisq(QIVW,DF)
Qind     = W1*(BIV - Bhat)^2

# IVWResults (Est., S.E, 95% CI, t-stat, p-value)

IVWResults = c(Bhat,SE,IVW_CI,Bhat/SE,IVW_p)

# IVWQ (Q, p-value,overdispersion)

IVWQ             = c(QIVW,Qp,phi_IVW)

return(list(IVWResults=IVWResults,IVWQ=IVWQ,Qind=Qind))

}


# 1st order IVW

W1           = 1/(seBetaYG^2/BetaXG^2)
IVW_1stOrder = IVWFunc(W=W1)

# 2nd Order IVW

W2           = 1/(seBetaYG^2/BetaXG^2 + (BetaYG^2)*seBetaXG^2/BetaXG^4)
IVW_2ndOrder = IVWFunc(W=W2)


# Iterative IVW

Diff  = 1
Bhat1 = 0 
print("Iterative IVW: Iteration #")
count = 0
while(Diff >= tol){
W3    = 1/(seBetaYG^2/BetaXG^2 + (Bhat1^2)*seBetaXG^2/BetaXG^2)
new   = IVWFunc(W=W3)
Diff  = abs(Bhat1 - new$IVWResults[1]) 
Bhat1 = new$IVWResults[1]
count = count+1
print(count)
}

IVW_Iterative = new 

lb = IVW_Iterative$IVWResults[1] - 5*IVW_Iterative$IVWResults[2]
ub = IVW_Iterative$IVWResults[1] + 5*IVW_Iterative$IVWResults[2]

#################
# Exact weights #
#################

PL2 = function(a){
b = a[1]
w = 1/((phi)*seBetaYG^2/BetaXG^2 + (b^2)*seBetaXG^2/BetaXG^2)
q =  sum(w*(BIV - b)^2)
}

PLfunc = function(a){
phi    = a[1]
PL2    = function(a){
beta   = a[1]
w      = 1/(phi*seBetaYG^2/BetaXG^2 + (beta^2)*seBetaXG^2/BetaXG^2)
q      =  (sum(w*(BIV - beta)^2))
}
b  = optimize(PL2,interval=c(lb,ub))$minimum 
w    = 1/(phi*seBetaYG^2/BetaXG^2 + (b^2)*seBetaXG^2/BetaXG^2)
q    =  (sum(w*(BIV - b)^2) - DF)^2
}

BootVar = function(sims=1000){
B = NULL ; pp=NULL
for(hh in 1:sims){
L      = length(BetaXG)
choice = sample(seq(1,L),L,replace=TRUE)
bxg    = BetaXG[choice] ; seX = seBetaXG[choice]
byg    = BetaYG[choice] ; seY = seBetaYG[choice]
BIV    = byg/bxg

W1        = 1/(seY^2/bxg^2)
BIVw1     = BIV*sqrt(W1)
sW1       = sqrt(W1)
IVWfitR1  = summary(lm(BIVw1 ~ -1+sW1))
phi_IVW1  = IVWfitR1$sigma^2
W2        = 1/(seY^2/bxg^2 + (byg^2)*seX^2/bxg^4)
BIVw2     = BIV*sqrt(W2)
sW2       = sqrt(W2)
IVWfitR2  = summary(lm(BIVw2 ~ -1+sW2))
phi_IVW2  = IVWfitR2$sigma^2

phi_IVW2 = max(1,phi_IVW2)
phi_IVW1 = max(1,phi_IVW1) 
lb       = IVWfitR1$coef[1] - 10*IVWfitR1$coef[2]
ub       = IVWfitR1$coef[1] + 10*IVWfitR1$coef[2]

PL2 = function(a){
b = a[1]
w = 1/((phi)*seY^2/bxg^2 + (b^2)*seX^2/bxg^2)
q =  sum(w*(BIV - b)^2)
}

PLfunc = function(a){
phi    = a[1]
PL2    = function(a){
beta   = a[1]
w      = 1/(phi*seY^2/bxg^2 + (beta^2)*seX^2/bxg^2)
q      =  (sum(w*(BIV - beta)^2))
}
b  = optimize(PL2,interval=c(-lb,ub))$minimum 
w    = 1/(phi*seY^2/bxg^2 + (b^2)*seX^2/bxg^2)
q    =  (sum(w*(BIV - b)^2) - DF)^2
}
phi    = optimize(PLfunc,interval=c(phi_IVW2,phi_IVW1+0.001))$minimum 
B[hh]  = optimize(PL2,interval=c(lb,ub))$minimum 
}
se   = sd(B)
mB   = mean(B)
return(list(mB=mB,se=se))
}

CIfunc = function(){
z = qt(df=DF, 0.975)
z2 = 2*(1-pnorm(z))

PL3 = function(a){
b = a[1]
w = 1/(seBetaYG^2/BetaXG^2 + (b^2)*seBetaXG^2/BetaXG^2)
q    =  (sum(w*(BIV - b)^2) - qchisq(1-z2,DF))^2
}

lb = Bhat - 10*SE
ub = Bhat + 10*SE

low   = optimize(PL3,interval=c(lb,Bhat))$minimum
high  = optimize(PL3,interval=c(Bhat,ub))$minimum
CI    = c(low,high)
return(list(CI=CI))
}


#######################
# Fixed effect model  #
# and Exact Q test    #
#######################

phi       = 1
Bhat      = optimize(PL2,interval=c(-2,2))$minimum 
W         = 1/(seBetaYG^2/BetaXG^2 + (Bhat^2)*seBetaXG^2/BetaXG^2)
SE        = sqrt(1/sum(W))
CI        = CIfunc()

# Qtest

QIVW      = sum(W*(BIV - Bhat)^2)
Qp        = 1-pchisq(QIVW,DF)
Qind      = W*(BIV - Bhat)^2
ExactQ    = c(QIVW,Qp)
ExactQind = Qind
# Estimation (fixed effects)
# point estimate, se, t-stat, p-value)

FE_EXACT    = c(Bhat,SE,CI$CI,Bhat/SE,2*(1-pt(abs(Bhat/SE),DF)))

# Estimation (random effects)

phi_IVW1 = IVW_1stOrder$IVWQ[3]
phi_IVW2 = IVW_2ndOrder$IVWQ[3]
phi_IVW2 = max(1,phi_IVW2)
phi_IVW1 = max(1,phi_IVW1)+0.001

Lb = Bhat - 10*SE
Ub = Bhat + 10*SE

phi       = optimize(PLfunc,interval=c(phi_IVW2,phi_IVW1))$minimum 
Bhat      = optimize(PL2,interval=c(Lb,Ub))$minimum 
Boot      = BootVar()
SE        = Boot$se 

# point estimate, se, t-stat, p-value)

IVW_CI    = Bhat + c(-1,1)*qt(df=DF, 0.975)*SE
RE_EXACT  = c(Bhat,SE,IVW_CI,Bhat/SE,2*(1-pt(abs(Bhat/SE),DF)))

RESULTS = matrix(nrow = 5,ncol=7)
RESULTS[1,] = c(IVW_1stOrder$IVWResults,IVW_1stOrder$IVWQ[3])
RESULTS[2,] = c(IVW_2ndOrder$IVWResults,IVW_2ndOrder$IVWQ[3])
RESULTS[3,] = c(IVW_Iterative$IVWResults,IVW_Iterative$IVWQ[3])
RESULTS[4,] = c(FE_EXACT,1)
RESULTS[5,] = c(RE_EXACT,phi)


RESULTS = data.frame(RESULTS)
colnames(RESULTS)=c("Est","SE","CI_L","CI_U","t-value","p-value","phi")
rownames(RESULTS)=c("1st order","2nd order","Iterative","Exact (FE)","Exact (RE)")

QStats     = matrix(nrow = 4,ncol=2)
QStats[1,] = IVW_1stOrder$IVWQ[1:2]
QStats[2,] = IVW_2ndOrder$IVWQ[1:2]
QStats[3,] = IVW_Iterative$IVWQ[1:2]
QStats[4,] = ExactQ

QStats = data.frame(QStats)
colnames(QStats) = c("Q statistic","p-value")
rownames(QStats)=c("1st order","2nd order","Iterative","Exact")

Qcontribution = matrix(nrow=(DF+1),ncol=4)
Qcontribution[,1] = IVW_1stOrder$Qind
Qcontribution[,2] = IVW_2ndOrder$Qind
Qcontribution[,3] = IVW_Iterative$Qind
Qcontribution[,4] = ExactQind

Qcontribution = data.frame(Qcontribution)
colnames(Qcontribution) = c("1st order","2nd Order","Iterative","Exact")


#plot(Qcontribution[,1],ylim=range(Qcontribution),xlab="SNP",
#ylab=" Cochran Q contribution",cex.lab=1.2,main="Q-contribution")
#points(Qcontribution[,2],col="black",pch=19)
#points(Qcontribution[,3],col="red",pch=19)
#points(Qcontribution[,4],col="green",pch=19)

#legend("topright",c("1st order weights", "2nd order weights", 
#       "Iterative","Exact"),
#       pch=c(1,19,19,19),col=c("black","black","red","green"),cex=1.2)

#lines(c(0,(DF+1)),rep(qchisq(0.95,1),2),lwd=2,lty=1)
#lines(c(0,(DF+1)),rep(qchisq(0.99,1),2),lwd=2,lty=2)
#lines(c(0,(DF+1)),rep(qchisq((1-0.05/(DF+1)),1),2),lwd=2,lty=3)

#legend("topleft",c("5%", "1%","(5/L)%"),lty=c(1,2,3))

############################################
# Apply weighted median with exact weights #
############################################

weighted.median <- function(betaIV.in, weights.in) {
betaIV.order = betaIV.in[order(betaIV.in)]
weights.order = weights.in[order(betaIV.in)]
weights.sum = cumsum(weights.order)-0.5*weights.order
weights.sum = weights.sum/sum(weights.order)
below = max(which(weights.sum<0.5))
weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
(0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
return(weighted.est) }

weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
med = NULL
for(i in 1:1000){
betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
betaIV.boot = betaYG.boot/betaXG.boot
med[i] = weighted.median(betaIV.boot, weights.in)
}
return(sd(med)) }

# Exact weight
W3 = 1/(phi*seBetaYG^2/BetaXG^2 + (Bhat^2)*seBetaXG^2/BetaXG^2)

betaIV   = BetaYG/BetaXG 
weights  = W3#(seBetaYG/BetaXG)^-2 
betaWM   = weighted.median(betaIV, weights) # weighted median estimate
sebetaWM = weighted.median.boot(BetaXG, BetaYG, seBetaXG, seBetaYG, weights) 
t     = betaWM/sebetaWM
p     = 2*(1-pt(abs(t),length(BetaYG)-1))
WMresults = data.frame(Estimate=betaWM,Std.Error=sebetaWM,t,p)
CI.WM     = betaWM + c(-1,1)*sebetaWM*qt(df=DF, 0.975)




return(list(RESULTS=RESULTS,QStats=QStats,
       Qcontribution=Qcontribution,meanF=mF,WM=WMresults))

}
