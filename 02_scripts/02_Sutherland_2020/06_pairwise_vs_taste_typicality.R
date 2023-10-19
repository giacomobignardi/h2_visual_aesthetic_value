#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description: repeat procedure in highlighted in 01_Germine_2015/06_pairwise_vs_taste_typicality.R
#Program: phenotyping_PCA ------------------------------------------------------------------------------------------------------------------------------

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
BioMetric  = read_csv(sprintf("%s/%s/02_Sutherland_2020/05_twin_multivariate_val.csv",wdOA,wdOA_output))
InterR_random = read_csv(sprintf("%s/%s/02_Sutherland_2020/03_pairwise_agreement_Sutherland2020.csv", wdOA,wdOA_output))
InterR_twin = InterR_random%>%filter(Pairs != "UR")

InterR_twin = InterR_random%>%
  filter(Pairs != "UR")%>%
  mutate(refCode = substr(refCode,1,nchar(refCode)-3))%>%
  rename(Sub = "refCode")

#calculate delta mm2 as the magnitude of the distance between two members of a pair
deltaMm2_twin = BioMetric%>%
  mutate(delta_mm2_SC = abs(SC_B - SC_A),
         delta_mm2_FA = abs(FA_B - FA_A)
  )%>%
  dplyr::select(Sub, delta_mm2_SC,delta_mm2_FA)

interr2deltamm2_twins = as_tibble(merge(InterR_twin,deltaMm2_twin, by = c("Sub"), all = T))


#check if the two metric are correlated
summary(lm(z_twin_PP~delta_mm2_FA,interr2deltamm2_twins[interr2deltamm2_twins$Domain == "faces val.",]))
lm(scale(z_twin_PP)~scale(delta_mm2_FA),interr2deltamm2_twins[interr2deltamm2_twins$Domain == "faces val.",])

summary(lm(z_twin_PP~delta_mm2_SC,interr2deltamm2_twins[interr2deltamm2_twins$Domain == "scenes val.",]))
lm(scale(z_twin_PP)~scale(delta_mm2_SC),interr2deltamm2_twins[interr2deltamm2_twins$Domain == "scenes val.",])

####VIZUALIZE####
p_dif_scenes_val = ggplot(interr2deltamm2_twins%>%filter(Domain == "scenes val."), aes(x = delta_mm2_SC, y= z_twin_PP )) + 
  geom_point(alpha =.5) +
  ylab("pa scenes val.") +
  xlab(expression(abs(Delta[paste(t-t,sep =" ",scenes,sep =" ",val)]))) +
  xlim(c(0,1.5))+
  ylim(c(0,2.2))+
  theme_classic(base_size = 16)
p_dif_faces_val =  ggplot(interr2deltamm2_twins%>%filter(Domain == "faces val."), aes(x = delta_mm2_FA, y= z_twin_PP )) + 
  geom_point(alpha =.5) +
  ylab("pa faces val.") +
  xlab(expression(abs(Delta[paste(t-t,sep =" ",faces,sep =" ",val)]))) +
  xlim(c(0,1.5))+
  ylim(c(0,2.2))+
  theme_classic(base_size = 16)

#SAVE###
save(p_dif_scenes_val,p_dif_faces_val, file = sprintf("%s/%s/02_Sutherland_2020/06_sup_FS4_tt2pairwise_val.rData",wdOA,wdOA_output))
   