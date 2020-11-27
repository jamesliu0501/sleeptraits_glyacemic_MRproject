# Overview
# This file contains the R script of the two-sample Mendelian randomization analyses (2SMRmain (IVW), 2SMRsensitivity1 (WM), 2SMRsensitivity2 (MR-Egger))
# Here we take sleep duration and insomnia vs HbA1c and as examples  (R codes of other sleep traits are similar)
# The id number of glucose (Dupuis J et al 2010, PMC: 3018764) from the MAGIC was 757 in the MRbase (version: "MRCIEU/TwoSampleMR@0.4.26")
# To be noticed, the effect allele frequency of the outcome summary data was misinterpreted in the “TwoSampleMR” package (version: "MRCIEU/TwoSampleMR@0.4.26")
# Therefore, we manually corrected the effect allele frequency of all the merged palindromic SNPs according to the data (i.e., minor allele frequency) downloaded from the MAGIC and the information of minor allele obtained from the ensembl of EUR population (https://www.ensembl.org/index.html)


# We used the old version of TwoSample MR 
# Because the updated version has moved the HbA1c GWAS from the MAGIC given the effect allele frequency was misinterpreted.
install.packages("devtools")
library(devtools)
devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")
devtools::install_github('MarkEdmondson1234/googleAuthR@v0.8.1')

library(TwoSampleMR)
library(dplyr)
library(MRInstruments)

rm(list = ls())
# Continuous exposure (e.g., sleep duration) 
# Sleep duration vs HbA1c
# Exposure summary statistics 
sleepduration <- read_exposure_data("sleepduration_exposure.csv", 
                                    sep=",", effect_allele = "A1", other_allele = "A2", eaf = "EAF", samplesize_col = "n" )

# To convert betaXG into hour unit (the original GWAS is in minute)
sleepduration$beta.exposure<-sleepduration$beta.exposure/60
sleepduration$se.exposure<-sleepduration$se.exposure/60

# Outcome 
# HbA1c
hba1c <- extract_outcome_data(sleepduration$SNP, outcome= "758") 
hba1c$effect_allele.outcome<-toupper(hba1c$effect_allele.outcome)
hba1c$other_allele.outcome<-toupper(hba1c$other_allele.outcome)

hba1c_outcome<-select(hba1c, 
                      SNP,
                      effect_allele.outcome, other_allele.outcome, 
                      beta.outcome, se.outcome, eaf.outcome)


# The MRbase misinterpreted the minor allele frequency (MAF) as the effect allele frequency (EAF)
# Therefore, we need to manually correct it as below
hba1c_outcome$maf.outcome<-hba1c_outcome$eaf.outcome 
hba1c_outcome$eaf.outcome<-NA


## To identify the palindromic SNP
sleepduration_dat <- harmonise_data(sleepduration, hba1c, action=3)
pSNP<-sleepduration_dat[sleepduration_dat$palindromic==TRUE,]
pSNP_only<-select(pSNP, SNP, 
                  proxy.outcome) # only identify the SNPs here. If we use this file to check the EAF, it would be wrong. Because the harominse function will flip the "EAF" based on the exposure file
pSNP_only$minor_allele.outcome<-NA
write.csv(pSNP_only, "sleepduration_hba1c_palindromic.csv",row.names = F)


## Check the minor allele manually from the ensembl of the EUR population (https://www.ensembl.org/index.html)
pSNP_updated<-read.csv("sleepduration_hba1c_palindromic_updated.csv")
pSNP_updated$palindromic=TRUE
hba1c_outcome_updated<-full_join(hba1c_outcome, pSNP_updated, by="SNP")


# Update eaf.outcome according to MAF and the minor allele
hba1c_outcome_updated$eaf.outcome<-ifelse(hba1c_outcome_updated$effect_allele.outcome == hba1c_outcome_updated$minor_allele.outcome,
                                          hba1c_outcome_updated$maf.outcome , 1 - hba1c_outcome_updated$maf.outcome)


# Update the outcome dataset and harominze again
outcomedata_updated<-select(hba1c_outcome_updated, 
                            SNP, 
                            effect_allele.outcome, other_allele.outcome, 
                            beta.outcome, se.outcome, eaf.outcome)

# Format the updated outcome dataset
outcomedata_updated<-format_data(outcomedata_updated, type = "outcome",
                                 snp_col="SNP",
                                 beta_col="beta.outcome", se_col="se.outcome", 
                                 effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome",
                                 eaf_col = "eaf.outcome")

write.csv(outcomedata_updated, "sleepduration_hba1c_outcome_updated.csv",row.names = F)

# Haromise again
final_haromized_data<-harmonise_data(sleepduration, outcomedata_updated, action=2)

write.csv(final_haromized_data, "sleepduration_hba1c_final_haromized_data.csv",row.names = F)



# Transform the units 
final_haromized_data<-read.csv("sleepduration_hba1c_final_haromized_data.csv")

# Update the outcome into SD unit (for HbA1c, 1SD = 0.5342848% as shown in the MRbase)
final_haromized_data$beta.outcome<-final_haromized_data$beta.outcome/0.5342848
final_haromized_data$se.outcome<-final_haromized_data$se.outcome/0.5342848

# Run MR regression and obtain the 2SMR estimates
sleepduration_hba1c_dat_mr <- mr(final_haromized_data) 
nsnp<-sleepduration_hba1c_dat_mr[1,6]

ivw_b<-sleepduration_hba1c_dat_mr[3,7]
ivw_lci<-ivw_b - 1.96*sleepduration_hba1c_dat_mr[3,8]
ivw_uci<-ivw_b + 1.96*sleepduration_hba1c_dat_mr[3,8]
ivw_p<-sleepduration_hba1c_dat_mr[3,9]

wm_b<-sleepduration_hba1c_dat_mr[2,7]
wm_lci<-wm_b - 1.96*sleepduration_hba1c_dat_mr[2,8]
wm_uci<-wm_b + 1.96*sleepduration_hba1c_dat_mr[2,8]
wm_p<-sleepduration_hba1c_dat_mr[2,9]

egger_b<-sleepduration_hba1c_dat_mr[1,7]
egger_lci<-egger_b - 1.96*sleepduration_hba1c_dat_mr[1,8]
egger_uci<-egger_b + 1.96*sleepduration_hba1c_dat_mr[1,8]
egger_p<-sleepduration_hba1c_dat_mr[1,9]


egger_intercept_p<-mr_pleiotropy_test(final_haromized_data)[1,7]















rm(list = ls())
# Binary exposure (e.g., insomnia) (we would use the BOLT-LMM method to convert logOR to beta)
# Insomnia vs HbA1c
# Exposure summary statistics 
insomnia <- read_exposure_data("insomnia_exposure.csv", 
                                    sep=",", effect_allele = "A1", other_allele = "A2", eaf = "EAF" , samplesize_col = "n")

# Outcome  
# Hba1c
hba1c <- extract_outcome_data(insomnia2$SNP, outcome= "758") 
hba1c$effect_allele.outcome<-toupper(hba1c$effect_allele.outcome)
hba1c$other_allele.outcome<-toupper(hba1c$other_allele.outcome)

hba1c_outcome<-select(hba1c, 
                      SNP,
                      effect_allele.outcome, other_allele.outcome, 
                      beta.outcome, se.outcome, eaf.outcome)


# The MRbase misinterpreted the minor allele frequency (MAF) as the effect allele frequency (EAF)
# Therefore, we need to manually correct it as below
hba1c_outcome$maf.outcome<-hba1c_outcome$eaf.outcome 
hba1c_outcome$eaf.outcome<-NA


## To identify the palindromic SNP
insomnia2_dat <- harmonise_data(insomnia2, hba1c, action=3)
pSNP<-insomnia2_dat[insomnia2_dat$palindromic==TRUE,]
pSNP_only<-select(pSNP, SNP, 
                  proxy.outcome) # only identify the SNPs here. If we use this file to check the EAF, it would be wrong. Since the harominse function will flip the "EAF" based on the exposure file
pSNP_only$minor_allele.outcome<-NA
write.csv(pSNP_only, "insomnia_hba1c_palindromic.csv",row.names = F)


## Check the minor allele manually from the ensembl of the EUR population (https://www.ensembl.org/index.html)
pSNP_updated<-read.csv("insomnia_hba1c_palindromic_updated.csv")
pSNP_updated$palindromic=TRUE
hba1c_outcome_updated<-full_join(hba1c_outcome, pSNP_updated, by="SNP")


# Update eaf.outcome according to MAF and the minor allele
hba1c_outcome_updated$eaf.outcome<-ifelse(hba1c_outcome_updated$effect_allele.outcome == hba1c_outcome_updated$minor_allele.outcome,
                                          hba1c_outcome_updated$maf.outcome , 1 - hba1c_outcome_updated$maf.outcome)


# Update the outcome dataset and harominze again
outcomedata_updated<-select(hba1c_outcome_updated, 
                            SNP, 
                            effect_allele.outcome, other_allele.outcome, 
                            beta.outcome, se.outcome, eaf.outcome)

# Format the updated outcome dataset
outcomedata_updated<-format_data(outcomedata_updated, type = "outcome",
                                 snp_col="SNP",
                                 beta_col="beta.outcome", se_col="se.outcome", 
                                 effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome",
                                 eaf_col = "eaf.outcome")

write.csv(outcomedata_updated, "insomnia_hba1c_outcome_updated.csv",row.names = F)



# Haromise again
final_haromized_data<-harmonise_data(insomnia2, outcomedata_updated, action=2)
write.csv(final_haromized_data, "insomnia_hba1c_final_haromized_data.csv",row.names = F)



# Transform the units 
final_haromized_data<-read.csv("insomnia_hba1c_final_haromized_data.csv")

# Update the outcome into SD unit (for HbA1c, 1SD = 0.5342848% as shown in the MRbase)
final_haromized_data$beta.outcome<-final_haromized_data$beta.outcome/0.5342848
final_haromized_data$se.outcome<-final_haromized_data$se.outcome/0.5342848


# Use the BOLT-LMM method to convert the results of the SNP-binary sleep trait from the multiplicative log odds scale to a difference in risk scale 
# β=log⁡OR* μ*(1-μ),se= 〖se〗_log⁡OR * μ*(1-μ), with μ=  n_case⁄((n_case+ n_control)).
# (https://data.bris.ac.uk/data/dataset/pnoat8cxo0u52p6ynfaekeigi)

insomnia_bolt<-read.csv("/insomnia2_exposure.csv")
u<-(insomnia_bolt$ncase.exposure / insomnia_bolt$n) [1]
final_haromized_data$beta.exposure<-final_haromized_data$beta.exposure * u * (1-u)
final_haromized_data$se.exposure<-final_haromized_data$se.exposure * u * (1-u)


# Run MR regression and obtain the 2SMR estimates
insomnia2_hba1c_dat_mr <- mr(final_haromized_data) 
nsnp<-insomnia2_hba1c_dat_mr[1,6]

ivw_b<-insomnia2_hba1c_dat_mr[3,7]
ivw_lci<-ivw_b - 1.96*insomnia2_hba1c_dat_mr[3,8]
ivw_uci<-ivw_b + 1.96*insomnia2_hba1c_dat_mr[3,8]
ivw_p<-insomnia2_hba1c_dat_mr[3,9]

wm_b<-insomnia2_hba1c_dat_mr[2,7]
wm_lci<-wm_b - 1.96*insomnia2_hba1c_dat_mr[2,8]
wm_uci<-wm_b + 1.96*insomnia2_hba1c_dat_mr[2,8]
wm_p<-insomnia2_hba1c_dat_mr[2,9]

egger_b<-insomnia2_hba1c_dat_mr[1,7]
egger_lci<-egger_b - 1.96*insomnia2_hba1c_dat_mr[1,8]
egger_uci<-egger_b + 1.96*insomnia2_hba1c_dat_mr[1,8]
egger_p<-insomnia2_hba1c_dat_mr[1,9]


egger_intercept_p<-mr_pleiotropy_test(final_haromized_data)[1,7]










