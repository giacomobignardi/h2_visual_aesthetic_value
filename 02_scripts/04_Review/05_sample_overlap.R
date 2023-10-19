#Author: Giacomo Bignardi
#Date: 19-09-2023
# Last modifed: 2023-10-10
#
#Description: 
#Estimate sample overlap between the Germine et al. and the Sutherland et al. study
#Program: S11_sample_overlap------------------------------------------------------------------------------------------------------------------------------

#Create a sample of 40000 individuals (Australian Twin Registry based on Hopper et al., 2013)
ATR_sim = c(1:40000)

#specify sample size of the Germine and the Sutherland studies as reported in main text
Germine_sim = 1115 + 432
Sutherland_sim = 815 + 416

#Draw two random sample from the ATR of the size of the Germine and the Sutherland studies
set.seed(42)
ovrlp_est = c()
for (i in 1:10000){
  sample1 = sample(ATR_sim,Germine_sim)
  sample2 = sample(ATR_sim,Sutherland_sim)
  ovrlp_est_i = sum(sample1 %in% sample2)
  ovrlp_est = c(ovrlp_est,ovrlp_est_i)
}

#mean overlap across the two samples
mean_overlapp = mean(ovrlp_est)
sd_overlapp = sd(ovrlp_est)
