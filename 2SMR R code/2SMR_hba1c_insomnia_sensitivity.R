# Overview
# This file contains the R script of the sensitivity analysis two-sample Mendelian randomization of HbA1c with insomnia 
# Here we take insomnia as an example  (R codes of other sleep traits are similar)

library(devtools)  
library(TwoSampleMR)
library(MRInstruments)
library(dplyr)

rm(list = ls())

# Exposure data of Hba1c is from MAGIC Soranzo N 2010  unit =per % 
exposure<-read.csv("hba1c_magic.csv", header  = TRUE)

# Format exposure data
exp_dat<-format_data(exposure,type = "exposure",
                     snp_col="SNP",beta_col="beta",
                     se_col="se",eaf_col="eaf", samplesize_col = "n",
                     effect_allele_col="effect_allele",other_allele_col="other_allele",
                     pval_col="pval")

# Clumping
exp_dat<-clump_data(exp_dat,clump_r2 = 0.001) #n=11


# Outcome data of insomnia (Jansen's GWAS) is downloaded from  https://ctg.cncr.nl/software/summary_statistics
# This file includes GWAS summary statistics of insomnia in UK Biobank (N=386,533 individuals)
# Due to restrictions on sharing the data, the results do not include 23andMe participants

outcome<- fread(file="Insomnia_sumstats_Jansenetal.txt"
                , fill=T, header  = TRUE)[,c("SNP", "A1", "A2", "MAF", "OR", "SE", "P", "N")] # MAF is minor allele frequency
colnames(outcome) <- c("SNP","effect_allele","other_allele","maf","or","se", "pval", "n")  # Cannot identify EAF, but there is no palindromic SNP
exp_snp<-select(exp_dat, SNP)
sub.df<-merge(exp_snp, outcome, by="SNP")
sub.df["phenotype"]<-"insomnia"

# Generate logodds
sub.df$beta<-log(sub.df$or)

# Format outcome data 
sub.df$effect_allele<-as.character(sub.df$effect_allele)
sub.df$other_allele<-as.character(sub.df$other_allele)
outcome_dat<-format_data(sub.df, type = "outcome",
                         snp_col="SNP",beta_col="beta",phenotype_col="phenotype",
                         se_col="se", samplesize_col = "n", pval_col="pval",
                         effect_allele_col="effect_allele",other_allele_col="other_allele")




# Harmonise data, 
#1 Assume all alleles are presented on the forward strand
#2 Try to infer the forward strand alleles using allele frequency information
#3 Correct the strand for non-palindromic SNPs, but drop all palindromic SNPs
dat<-harmonise_data(exposure_dat = exp_dat,
                    outcome_dat = outcome_dat,
                    action = 1) # Cannot identify EAF, but there is no palindromic SNP


# Run MR regressions 
hba1c_insomnia_mr <- mr(dat) 
nsnp<-hba1c_insomnia_mr[1,6]


ivw_b<-exp(hba1c_insomnia_mr[3,7])
ivw_lci<-exp(hba1c_insomnia_mr[3,7] - 1.96*hba1c_insomnia_mr[3,8])
ivw_uci<-exp(hba1c_insomnia_mr[3,7] + 1.96*hba1c_insomnia_mr[3,8])
ivw_p<-hba1c_insomnia_mr[3,9]

wm_b<-exp(hba1c_insomnia_mr[2,7])
wm_lci<-exp(hba1c_insomnia_mr[2,7] - 1.96*hba1c_insomnia_mr[2,8])
wm_uci<-exp(hba1c_insomnia_mr[2,7] + 1.96*hba1c_insomnia_mr[2,8])
wm_p<-hba1c_insomnia_mr[2,9]

egger_b<-exp(hba1c_insomnia_mr[1,7])
egger_lci<-exp(hba1c_insomnia_mr[1,7] - 1.96*hba1c_insomnia_mr[1,8])
egger_uci<-exp(hba1c_insomnia_mr[1,7] + 1.96*hba1c_insomnia_mr[1,8])
egger_p<-hba1c_insomnia_mr[1,9]


egger_intercept_p<-mr_pleiotropy_test(dat)[1,7]

