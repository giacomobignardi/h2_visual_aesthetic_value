#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
# First create a list of exclusion criteria then create a set of dataframe to be used in following analysis
# intra-rater -> exclusion criteria based on Vessel et al. (intrarate <.5)
# Biometric -> to be used for following biometric modeling
# MLM -> to be used to fit the multilevel model to do VCA
# MM@ -> to be used to compute mean minus 2 scores

#Program: twinDfs ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)
library(umx)
library(patchwork)


#clean working enviroment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"
#load dataFrames:
BioMetric  = read_csv(sprintf("%s/%s/01_Germine_2015/05_twin_univariate.csv",wdOA,wdOA_output))
InterR_random = read_csv(sprintf("%s/%s/01_Germine_2015/03_pairwise_agreement_Germine2015.csv", wdOA,wdOA_output))
InterR_twin = InterR_random%>%filter(Pairs != "UR")

#calculate delta mm2 as the magnitude of the distance between two members of a pair
deltaMm2_twin = BioMetric%>%
  mutate(delta_mm2_z = abs(mm2_z_2 - mm2_z_1),
         Domain= fct_recode(category, faces = "FA_TOT", scenes = "SC", abstract = "AO"))%>%
  dplyr::select(Domain,FamId,delta_mm2_z)

interr2deltamm2_twins = as_tibble(merge(InterR_twin,deltaMm2_twin, by = c("FamId", "Domain"), all = T))

#check if the two metric are correlated
summary(lm(z_twin_PP~delta_mm2_z,interr2deltamm2_twins[interr2deltamm2_twins$Domain == "faces",]))
lm(scale(delta_mm2_z)~scale(z_twin_PP),interr2deltamm2_twins[interr2deltamm2_twins$Domain == "faces",])

summary(lm(z_twin_PP~delta_mm2_z,interr2deltamm2_twins[interr2deltamm2_twins$Domain == "scenes",]))
lm(scale(delta_mm2_z)~scale(z_twin_PP),interr2deltamm2_twins[interr2deltamm2_twins$Domain == "scenes",])

summary(lm(z_twin_PP~delta_mm2_z,interr2deltamm2_twins[interr2deltamm2_twins$Domain == "abstract",]))
lm(scale(delta_mm2_z)~scale(z_twin_PP),interr2deltamm2_twins[interr2deltamm2_twins$Domain == "abstract",])

####VIZUALIZE####
p_dif_scenes = ggplot(interr2deltamm2_twins%>%filter(Domain == "scenes"), aes(x = delta_mm2_z, y= z_twin_PP )) + 
  geom_point(alpha =.5) +
  ylab("pa scenes") +
  xlab(expression(abs(Delta[paste(t-t,sep =" ",scenes)]))) +
  xlim(c(0,1.5))+
  ylim(c(-0.5,2.2))+
  theme_classic(base_size = 16)
p_dif_faces =  ggplot(interr2deltamm2_twins%>%filter(Domain == "faces"), aes(x = delta_mm2_z, y= z_twin_PP )) + 
  geom_point(alpha =.5) +
  ylab("pa faces") +
  xlab(expression(abs(Delta[paste(t-t,sep =" ",faces)]))) +
  xlim(c(0,1.5))+
  ylim(c(-0.5,2.2))+
  theme_classic(base_size = 16)
p_dif_abstract = ggplot(interr2deltamm2_twins%>%filter(Domain == "abstract"), aes(x = delta_mm2_z, y= z_twin_PP )) + 
  geom_point(alpha =.5) +
  ylab("pa abstract") +
  xlab(expression(abs(Delta[paste(t-t,sep =" ",abstract)]))) +
  xlim(c(0,1.5))+
  ylim(c(-0.5,2.2))+
  theme_classic(base_size = 16)

save(p_dif_scenes,p_dif_faces,p_dif_abstract, file = sprintf("%s/%s/01_Germine_2015/06_sup_FS4_tt2pairwise.rData",wdOA,wdOA_output))

