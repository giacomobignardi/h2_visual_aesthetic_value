#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
#1. Create and prepare Evaluation-bias and Taste-typicality score for later analysis
#2. Carry out Principal Component Analysis on aesthetic ratings
#3. Correlate PCs with Evaluation-bias and Taste-typicality scores
#4. Prepare df for univariate and multivariate twin modeling
#Program: 05_phenotyping_PCA ------------------------------------------------------------------------------------------------------------------------------

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
BioMetric  = read_csv(sprintf("%s/%s/01_Germine_2015/01_CTD_Germine2015.csv", wdOA,wdOA_output))
Fac_MM2_Output = read_csv(sprintf("%s/%s/01_Germine_2015/04_taste_typicality_faces_Germine2015.csv", wdOA,wdOA_output))
Sce_MM2_Output = read_csv(sprintf("%s/%s/01_Germine_2015/04_taste_typicality_scenes_Germine2015.csv", wdOA,wdOA_output))
Obj_MM2_Output = read_csv(sprintf("%s/%s/01_Germine_2015/04_taste_typicality_abstracts_Germine2015.csv", wdOA,wdOA_output))


#EVALUATION-BIAS#####
#apply exclusion criteria
#average repeated measure first
Fac_BioMetric_avg = BioMetric%>%
  filter(category == "FA_TOT", intraR_FA >=.5)%>%
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Fac_BioMetric = Fac_BioMetric_avg %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Sce_BioMetric_avg  = BioMetric%>%
  filter(category == "SC", intraR_SC >=.5)%>%
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Sce_BioMetric = Sce_BioMetric_avg %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Obj_BioMetric_avg = BioMetric%>%
  filter(category == "AO", intraR_AO >=.5)%>%
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Obj_BioMetric = Obj_BioMetric_avg %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

#MM2####
Fac_MM2_Output = Fac_MM2_Output%>%
  mutate(category = "FA_TOT")
Sce_MM2_Output = Sce_MM2_Output%>%
  mutate(category = "SC")
Obj_MM2_Output = Obj_MM2_Output%>%
  mutate(category = "AO")

#bind al together (both overall pleasantess and taste typicality)
BioMetricThin = rbind(Fac_BioMetric,Obj_BioMetric,Sce_BioMetric)
MM2_Output = rbind(Fac_MM2_Output,Sce_MM2_Output,Obj_MM2_Output)

#####Outliers####
##PCA is  sensitive to extreme outliers: as such remove extreme outliers for specific domain from subsequent analysis
#exclude extreme outliers
Outliers_mm2 = MM2_Output%>%
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub)))   %>% 
  group_by(category, SibId)%>%
  rstatix::identify_outliers(mm2_z)%>%
  filter(is.extreme == T) %>% 
  select(category, Sub) 
Outliers_mm2_FA = Outliers_mm2 %>% filter(category == "FA_TOT") %>% pull(Sub)#2
Outliers_mm2_SC = Outliers_mm2 %>% filter(category == "SC") %>% pull(Sub)#2
Outliers_p = BioMetricThin%>%
  group_by(category, SibId)%>%
  rstatix::identify_outliers(avgP)%>%
  filter(is.extreme == T) %>% 
  filter(category != "DO") %>% 
  pull(Sub)

#Vizualize distributions####
#Overall pleasantess
p_avgP_faces = ggplot(BioMetricThin%>%
                        filter(category=="FA_TOT")) + 
  geom_histogram(bins = 50, aes(x = avgP), position="identity", colour='white') + 
  labs(x = "o-p faces")+
  scale_y_continuous(limits= c(0,110),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(2,6))+
  theme_classic(base_size=16)

p_avgP_scenes =ggplot(BioMetricThin%>%
                        filter(category=="SC")) + 
  geom_histogram(bins = 50, aes(x = avgP), position="identity", colour='white') + 
  labs(x = "o-p scenes")+
  scale_y_continuous(limits= c(0,110),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(2,6))+
  theme_classic(base_size=16)

p_avgP_abstracts =ggplot(BioMetricThin%>%
                           filter(category=="AO")) + 
  geom_histogram(bins = 50, aes(x = avgP), position="identity", colour='white') + 
  labs(x = "o-p abstracts")+
  scale_y_continuous(limits= c(0,110),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(2,6))+
  theme_classic(base_size=16)

p_avgP_faces| p_avgP_scenes| p_avgP_abstracts

#mm2
p_mm2_faces = ggplot(MM2_Output%>%
                       filter(category=="FA_TOT") %>% 
                       filter(!( Sub  %in% c(Outliers_mm2_FA)))) + 
  geom_histogram(bins = 50, aes(x = mm2_z), position="identity", colour='white') + 
  labs(x = "t-t faces")+
  scale_x_continuous(limits= c(-.6,2.5),breaks = c(0,1,2))+
  scale_y_continuous(limits= c(0,280),breaks = c(0,80,160,240))+
  theme_classic(base_size=16)

p_mm2_scenes = ggplot(MM2_Output%>%
                        filter(category=="SC") %>% 
                        filter(!( Sub  %in% c(Outliers_mm2_SC)))) + 
  geom_histogram(bins = 50, aes(x = mm2_z), position="identity", colour='white') + 
  labs(x = "t-t scenes")+
  scale_x_continuous(limits= c(-.6,2.5),breaks = c(0,1,2))+
  scale_y_continuous(limits= c(0,240),breaks = c(0,80,160,240))+
  theme_classic(base_size=16)

p_mm2_abstracts = ggplot(MM2_Output%>%filter(category=="AO")) + 
  geom_histogram(bins = 50, aes(x = mm2_z), position="identity", colour='white') + 
  labs(x = "t-t asbtracts",
       subtitle = "c")+
  scale_y_continuous(limits= c(0,240),breaks = c(0,40,80,120))+
  scale_x_continuous(limits= c(-.6,2.5),breaks = c(0,1,2))+
  theme_classic(base_size=16)

p_mm2_faces|p_mm2_scenes|p_mm2_abstracts


#PCA####
set.seed(123)
#####compute PC####
Fac_BioMetric_pc_t1 = Fac_BioMetric_avg %>% 
  filter(SibId == 1) %>% 
  filter(!( FamId_2  %in% c(Outliers_mm2_FA))) %>% 
  select(FamId_2, Value ,Item_2) %>% 
  pivot_wider(values_from ="Value", names_from = "Item_2") #prepare wide dataframe (row individuals, column images)
psych::KMO(Fac_BioMetric_pc_t1[,-1])
Fac_PC_t1 = prcomp(Fac_BioMetric_pc_t1[,-1]) #Principal Componet Analysis
Fac_PCscore_t1 = cbind(Fac_BioMetric_pc_t1[,1],Fac_PC_t1$x) #Extract score per individual and bind with individual ID

Sce_BioMetric_pc_t1 = Sce_BioMetric_avg %>% 
  filter(SibId == 1) %>% 
  filter(!( FamId_2  %in% c(Outliers_mm2_SC))) %>% 
  select(FamId_2, Value ,Item_2) %>% 
  pivot_wider(values_from ="Value", names_from = "Item_2")
psych::KMO(Sce_BioMetric_pc_t1[,-1])
Sce_PC_t1 = prcomp(Sce_BioMetric_pc_t1[,-1])
Sce_PCscore_t1 = cbind(Sce_BioMetric_pc_t1[,1],Sce_PC_t1$x)

Obj_BioMetric_pc_t1 = Obj_BioMetric_avg %>% 
  filter(SibId == 1) %>% 
  select(FamId_2, Value ,Item_2) %>% 
  pivot_wider(values_from ="Value", names_from = "Item_2")
psych::KMO(Obj_BioMetric_pc_t1[,-1])
Obj_PC_t1 = prcomp(Obj_BioMetric_pc_t1[,-1])
Obj_PCscore_t1 = cbind(Obj_BioMetric_pc_t1[,1],Obj_PC_t1$x)

Fac_BioMetric_pc_t2 = Fac_BioMetric_avg %>% 
  filter(SibId == 2) %>% 
  filter(!( FamId_2  %in% c(Outliers_mm2_FA))) %>% 
  select(FamId_2, Value ,Item_2) %>% 
  pivot_wider(values_from ="Value", names_from = "Item_2")
psych::KMO(Fac_BioMetric_pc_t2[,-1])
Fac_PC_t2 = prcomp(Fac_BioMetric_pc_t2[,-1])
Fac_PCscore_t2 = cbind(Fac_BioMetric_pc_t2[,1],Fac_PC_t2$x)

Sce_BioMetric_pc_t2 = Sce_BioMetric_avg %>% 
  filter(SibId == 2) %>% 
  filter(!( FamId_2  %in% c(Outliers_mm2_SC))) %>% 
  select(FamId_2, Value ,Item_2) %>% 
  pivot_wider(values_from ="Value", names_from = "Item_2")
psych::KMO(Sce_BioMetric_pc_t2[,-1])
Sce_PC_t2 = prcomp(Sce_BioMetric_pc_t2[,-1])
Sce_PCscore_t2 = cbind(Sce_BioMetric_pc_t2[,1],Sce_PC_t2$x)

Obj_BioMetric_pc_t2 = Obj_BioMetric_avg %>% 
  filter(SibId == 2) %>% 
  select(FamId_2, Value ,Item_2) %>% 
  pivot_wider(values_from ="Value", names_from = "Item_2")
psych::KMO(Obj_BioMetric_pc_t2[,-1])
Obj_PC_t2 = prcomp(Obj_BioMetric_pc_t2[,-1])
Obj_PCscore_t2 = cbind(Obj_BioMetric_pc_t2[,1],Obj_PC_t2$x)

#####Proportion of variance####
#Variance explained
FA_PC_t1  = Fac_PC_t1$sdev^2 / sum(Fac_PC_t1$sdev^2)
SC_PC_t1  = Sce_PC_t1$sdev^2 / sum(Sce_PC_t1$sdev^2)
AO_PC_t1  = Obj_PC_t1$sdev^2 / sum(Obj_PC_t1$sdev^2)

FA_PC_t2  = Fac_PC_t2$sdev^2 / sum(Fac_PC_t2$sdev^2)
SC_PC_t2  = Sce_PC_t2$sdev^2 / sum(Sce_PC_t2$sdev^2)
AO_PC_t2  = Obj_PC_t2$sdev^2 / sum(Obj_PC_t2$sdev^2)

#total variance
PC12_FA_var_explained_t1_Germine = sum(FA_PC_t1[1:2])
PC12_SC_var_explained_t1_Germine = sum(SC_PC_t1[1:2])
PC12_AO_var_explained_t1_Germine =  sum(AO_PC_t1[1:2])

PC12_FA_var_explained_t2_Germine = sum(FA_PC_t2[1:2])
PC12_SC_var_explained_t2_Germine = sum(SC_PC_t2[1:2])
PC12_AO_var_explained_t2_Germine = sum(AO_PC_t2[1:2])

#####Correlation with facets####
#bind twin one and twin two and prepare to add info to dataframe
Fac_PCscore = rbind(Fac_PCscore_t1[,1:3], Fac_PCscore_t2[,1:3])
Sce_PCscore = rbind(Sce_PCscore_t1[,1:3], Sce_PCscore_t2[,1:3])
Obj_PCscore = rbind(Obj_PCscore_t1[,1:3], Obj_PCscore_t2[,1:3])

Fac_BioMetric_PC_avg = merge(Fac_BioMetric %>% filter(!( Sub  %in% c(Outliers_mm2_FA))),Fac_PCscore %>% rename(Sub   = "FamId_2"), by = c("Sub"), all = T)
Sce_BioMetric_PC_avg = merge(Sce_BioMetric %>% filter(!( Sub  %in% c(Outliers_mm2_SC))),Sce_PCscore %>% rename(Sub   = "FamId_2"), by = c("Sub"), all = T)
Obj_BioMetric_PC_avg = merge(Obj_BioMetric ,Obj_PCscore %>% rename(Sub   = "FamId_2"), by = c("Sub"), all = T)

Fac_BioMetric_PC_mm2 = merge(Fac_MM2_Output %>% filter(!( Sub  %in% c(Outliers_mm2_FA))),Fac_PCscore %>% rename(Sub   = "FamId_2"), by = c("Sub"), all = T);Fac_BioMetric_PC_mm2 =  Fac_BioMetric_PC_mm2 %>% mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub)))
Sce_BioMetric_PC_mm2 = merge(Sce_MM2_Output %>% filter(!( Sub  %in% c(Outliers_mm2_SC))),Sce_PCscore %>% rename(Sub   = "FamId_2"), by = c("Sub"), all = T);Sce_BioMetric_PC_mm2 =  Sce_BioMetric_PC_mm2 %>% mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub)))
Obj_BioMetric_PC_mm2 = merge(Obj_MM2_Output,Obj_PCscore %>% rename(Sub   = "FamId_2"), by = c("Sub"), all = T);Obj_BioMetric_PC_mm2 =  Obj_BioMetric_PC_mm2 %>% mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub)))


#Correlation between PCs and Facets of aesthetic value
PC_FA_cor_eb_t1_Germine = cor.test(-(Fac_BioMetric_PC_avg %>% filter(SibId == 1) %>% pull(PC1)),  Fac_BioMetric_PC_avg %>% filter(SibId == 1) %>% pull(avgP))
PC_SC_cor_eb_t1_Germine = cor.test(-(Sce_BioMetric_PC_avg %>% filter(SibId == 1) %>% pull(PC1)),  Sce_BioMetric_PC_avg %>% filter(SibId == 1) %>% pull(avgP))
PC_AO_cor_eb_t1_Germine = cor.test(-(Obj_BioMetric_PC_avg %>% filter(SibId == 1) %>% pull(PC1)),  Obj_BioMetric_PC_avg %>% filter(SibId == 1) %>% pull(avgP))
PC_FA_cor_eb_t2_Germine = cor.test(-(Fac_BioMetric_PC_avg %>% filter(SibId == 2) %>% pull(PC1)),  Fac_BioMetric_PC_avg %>% filter(SibId == 2) %>% pull(avgP))
PC_SC_cor_eb_t2_Germine = cor.test(-(Sce_BioMetric_PC_avg %>% filter(SibId == 2) %>% pull(PC1)),  Sce_BioMetric_PC_avg %>% filter(SibId == 2) %>% pull(avgP))
PC_AO_cor_eb_t2_Germine = cor.test(-(Obj_BioMetric_PC_avg %>% filter(SibId == 2) %>% pull(PC1)),  Obj_BioMetric_PC_avg %>% filter(SibId == 2) %>% pull(avgP))

PC_FA_cor_tt_t1_Germine = cor.test(-(Fac_BioMetric_PC_mm2 %>% filter(SibId == 1) %>% pull(PC2)),  Fac_BioMetric_PC_mm2 %>% filter(SibId == 1) %>% pull(mm2_z))
PC_SC_cor_tt_t1_Germine = cor.test(-(Sce_BioMetric_PC_mm2 %>% filter(SibId == 1) %>% pull(PC2)),  Sce_BioMetric_PC_mm2 %>% filter(SibId == 1) %>% pull(mm2_z))
PC_AO_cor_tt_t1_Germine = cor.test(-(Obj_BioMetric_PC_mm2 %>% filter(SibId == 1) %>% pull(PC2)),  Obj_BioMetric_PC_mm2 %>% filter(SibId == 1) %>% pull(mm2_z))
PC_FA_cor_tt_t2_Germine = cor.test((Fac_BioMetric_PC_mm2 %>% filter(SibId == 2) %>% pull(PC2)),  Fac_BioMetric_PC_mm2 %>% filter(SibId == 2) %>% pull(mm2_z)) #flipped
PC_SC_cor_tt_t2_Germine = cor.test((Sce_BioMetric_PC_mm2 %>% filter(SibId == 2) %>% pull(PC2)),  Sce_BioMetric_PC_mm2 %>% filter(SibId == 2) %>% pull(mm2_z)) #flipped
PC_AO_cor_tt_t2_Germine = cor.test((Obj_BioMetric_PC_mm2 %>% filter(SibId == 2) %>% pull(PC2)),  Obj_BioMetric_PC_mm2 %>% filter(SibId == 2) %>% pull(mm2_z)) #flipped

#REPORT:PCA####
save(PC_FA_cor_eb_t1_Germine, PC_SC_cor_eb_t1_Germine, PC_AO_cor_eb_t1_Germine,
     PC_FA_cor_eb_t2_Germine, PC_SC_cor_eb_t2_Germine, PC_SC_cor_eb_t2_Germine,
     PC_FA_cor_tt_t1_Germine, PC_SC_cor_tt_t1_Germine, PC_AO_cor_tt_t1_Germine,
     PC_FA_cor_tt_t2_Germine, PC_SC_cor_tt_t2_Germine, PC_AO_cor_tt_t2_Germine,
     PC12_FA_var_explained_t1_Germine,
     PC12_SC_var_explained_t1_Germine,
     PC12_AO_var_explained_t1_Germine,
     PC12_FA_var_explained_t2_Germine,
     PC12_SC_var_explained_t2_Germine,
     PC12_AO_var_explained_t2_Germine,
     file = sprintf("%s/%s/01_Germine_2015/05_report_PCA.Rdata",wdOA,wdOA_output))



#####F1e####
p_Fac_BioMetric_PC_avg_t1_Ger2015 = Fac_BioMetric_PC_avg %>% 
  filter(SibId == 1) %>% 
  ggplot(aes(y = -PC1,avgP)) + 
  xlim(1,7)+
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
       y = "(-)",
       x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Sce_BioMetric_PC_avg_t1_Ger2015 = Sce_BioMetric_PC_avg %>% 
  filter(SibId == 1) %>% 
  ggplot(aes(y =-PC1,avgP)) + 
  xlim(1,7)+
  geom_point(fill = "black",colour = "white",pch=21,  size = 1.5) + 
  labs(#title = "Scenes",
       y = "(-)",
       x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Obj_BioMetric_PC_avg_t1_Ger2015 = Obj_BioMetric_PC_avg %>% 
  filter(SibId == 1) %>% 
  ggplot(aes(y = -PC1,avgP)) + 
  xlim(1,7)+
  geom_point(fill = "black",colour = "white",pch=21,  size = 1.5) + 
  labs(#title = "Abstracts",
       y = "(-) PC1 score",
       x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Fac_BioMetric_PC_mm2_t1_Ger2015 = Fac_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibId == 1) %>% 
  ggplot(aes(y =-PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
       y = "(-)",
       x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Sce_BioMetric_PC_mm2_t1_Ger2015 = Sce_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibId == 1) %>% 
  ggplot(aes(y =-PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Scenes",
       y = "(-)",
       x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Obj_BioMetric_PC_mm2_t1_Ger2015 = Obj_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibId == 1) %>% 
  ggplot(aes(y =-PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21 , size = 1.5) + 
  labs(#title = "Abstracts",
       y = "(-)PC2 score",
       x = "taste-typicality") +
  theme_classic(base_size = 16) 

#Twin 2
p_Fac_BioMetric_PC_avg_t2_Ger2015 = Fac_BioMetric_PC_avg %>% 
  filter(SibId == 2) %>% 
  ggplot(aes(y = -PC1,avgP)) + 
  xlim(1,7)+
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
    y = "(-)",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Sce_BioMetric_PC_avg_t2_Ger2015 = Sce_BioMetric_PC_avg %>% 
  filter(SibId == 2) %>% 
  ggplot(aes(y =-PC1,avgP)) + 
  xlim(1,7)+
  geom_point(fill = "black",colour = "white",pch=21,  size = 1.5) + 
  labs(#title = "Scenes",
    y = "(-)",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Obj_BioMetric_PC_avg_t2_Ger2015 = Obj_BioMetric_PC_avg %>% 
  filter(SibId == 2) %>% 
  ggplot(aes(y = -PC1,avgP)) + 
  xlim(1,7)+
  geom_point(fill = "black",colour = "white",pch=21,  size = 1.5) + 
  labs(#title = "Abstracts",
    y = "(-)PC1 score",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Fac_BioMetric_PC_mm2_t2_Ger2015 = Fac_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibId == 2) %>% 
  ggplot(aes(y =PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
    y = "",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Sce_BioMetric_PC_mm2_t2_Ger2015 = Sce_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibId == 2) %>% 
  ggplot(aes(y =PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Scenes",
    y = "",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Obj_BioMetric_PC_mm2_t2_Ger2015 = Obj_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibId == 2) %>% 
  ggplot(aes(y =PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21 , size = 1.5) + 
  labs(#title = "Abstracts",
    y = "PC2 score",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 


#Save PC
save(p_Fac_BioMetric_PC_avg_t1_Ger2015, p_Fac_BioMetric_PC_avg_t2_Ger2015, p_Sce_BioMetric_PC_avg_t1_Ger2015, p_Sce_BioMetric_PC_avg_t2_Ger2015,p_Obj_BioMetric_PC_avg_t1_Ger2015, p_Obj_BioMetric_PC_avg_t2_Ger2015,
     p_Fac_BioMetric_PC_mm2_t1_Ger2015, p_Fac_BioMetric_PC_mm2_t2_Ger2015, p_Sce_BioMetric_PC_mm2_t1_Ger2015, p_Sce_BioMetric_PC_mm2_t2_Ger2015,p_Obj_BioMetric_PC_mm2_t1_Ger2015,p_Obj_BioMetric_PC_mm2_t2_Ger2015,
     file = sprintf("%s/%s/01_Germine_2015/05_sup_S5c_PC_Germine2015.Rdata",wdOA,wdOA_output))


####F1c & FS5####
#Plot PC1 PC2 score for abstract and overlay to pleasantess and taste-typicality scores
p_PC1_FA = 
  ggplot(Fac_BioMetric_PC_avg %>% filter(SibId == 1), aes(-PC1, -PC2, color = avgP)) +
  geom_point() +
  labs(x = paste0("(-)PC1 (",round(FA_PC_t1[1],2),"%)"),
       y = paste0("(-)PC2 (",round(FA_PC_t1[2],2),"%)"),
       color = "Evaluation-bias",
       subtitle = "Faces")+
  scale_color_viridis_c(option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC2_FA =ggplot(Fac_BioMetric_PC_mm2 %>% filter(SibId == 1), aes(-PC1, -PC2, color = mm2_z)) +
  geom_point() +
  labs(x = paste0("(-)PC1 (",round(FA_PC_t1[1],2),"%)"),
       y =  paste0("(-)PC2 (",round(FA_PC_t1[2],2),"%)"),
       color = "Taste-typicality",
       subtitle = "Faces")+
  scale_color_viridis_c(option = "cividis") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC1_SC = 
  ggplot(Sce_BioMetric_PC_avg %>% filter(SibId == 1), aes(-PC1, -PC2, color = avgP)) +
  geom_point() +
  labs(x = paste0("(-)PC1 (",round(SC_PC_t1[1],2),"%)"),
       y = paste0("(-)PC2 (",round(SC_PC_t1[2],2),"%)"),
       color = "Evaluation-bias",
       subtitle = "Scenes")+
  scale_color_viridis_c(option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC2_SC =ggplot(Sce_BioMetric_PC_mm2 %>% filter(SibId == 1), aes(-PC1, -PC2, color = mm2_z)) +
  geom_point() +
  labs(x = paste0("(-)PC1 (",round(SC_PC_t1[1],2),"%)"),
       y =  paste0("(-)PC2 (",round(SC_PC_t1[2],2),"%)"),
       color = "Taste-typicality",
       subtitle = "Scenes")+
  scale_color_viridis_c(option = "cividis") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC1_AO = 
  ggplot(Obj_BioMetric_PC_avg %>% filter(SibId == 1), aes(-PC1, -PC2, color = avgP)) +
  geom_point() +
  labs(x = paste0("(-)PC1 (",round(AO_PC_t1[1],2),"%)"),
       y = paste0("(-)PC2 (",round(AO_PC_t1[2],2),"%)"),
       color = "Evaluation-bias",
       subtitle = "Abstracts")+
  scale_color_viridis_c(option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC2_AO =ggplot(Obj_BioMetric_PC_mm2 %>% filter(SibId == 1), aes(-PC1, -PC2, color = mm2_z)) +
  geom_point() +
  labs(x = paste0("(-)PC1 (",round(AO_PC_t1[1],2),"%)"),
       y = paste0("(-)PC2 (",round(AO_PC_t1[2],2),"%)"),
       color = "Taste-typicality",
       subtitle = "Abstracts")+
  scale_color_viridis_c(option = "cividis") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

pdf(sprintf("%s/%s/05_F1c_PC_AO.pdf",wdOA,wdOA_ImageOutput),
    width = 9,
    height = 4)
(p_PC2_AO | p_PC1_AO) + plot_layout(guides = "collect")
dev.off()

#####FS5####
#scree plot
#take variance explained from each component (select first 10)
PC_t1_sp = rbind(
  FA_PC_t1_sp = data.frame(prop_variance = FA_PC_t1[1:10], component = rep(1:10), visual_domain = "faces"),
  SC_PC_t1_sp = data.frame(prop_variance = SC_PC_t1[1:10],component =  rep(1:10),  visual_domain = "scenes"),
  AO_PC_t1_sp =  data.frame(prop_variance = AO_PC_t1[1:10], component = rep(1:10),  visual_domain = "abstracts")
)

PC_t1_sp = PC_t1_sp %>% mutate(visual_domain = factor(visual_domain, levels= c("faces","scenes","abstracts")))
#take variance explained from each component (select first 10) for twin 2
PC_t2_sp = rbind(
  FA_PC_t2_sp = data.frame(prop_variance = FA_PC_t2[1:10], component = rep(1:10), visual_domain = "faces"),
  SC_PC_t2_sp = data.frame(prop_variance = SC_PC_t2[1:10],component =  rep(1:10),  visual_domain = "scenes"),
  AO_PC_t2_sp =  data.frame(prop_variance = AO_PC_t2[1:10], component = rep(1:10),  visual_domain = "abstracts")
)

PC_t2_sp = PC_t2_sp %>% mutate(visual_domain = factor(visual_domain, levels= c("faces","scenes","abstracts")))

#Scree plot
p_sp_1 = ggplot(data=PC_t1_sp, aes(x=component, y=prop_variance, color = visual_domain)) +
  geom_line()+
  scale_x_continuous(breaks=seq(1, 10,1))+
  ylim(c(0,0.5))+
  scale_color_brewer(palette = "Dark2")+
  geom_point() +
  theme_classic()
p_sp_2 = ggplot(data=PC_t2_sp, aes(x=component, y=prop_variance, color = visual_domain)) +
  geom_line()+
  scale_x_continuous(breaks=seq(1, 10,1))+
  ylim(c(0,0.5))+
  scale_color_brewer(palette = "Dark2")+
  geom_point() +
  theme_classic()

save(p_PC2_AO,
     p_PC1_AO,
     p_PC2_SC,
     p_PC1_SC,
     p_PC2_FA, 
     p_PC1_FA,
     p_sp_1,
     p_sp_2,
     file = sprintf("%s/%s/01_Germine_2015/05_sup_FS5_screeplot_PCA.Rdata",wdOA,wdOA_output))

#UNIVARIATE CTD####
#merge data in one
MM2_Output[MM2_Output$Sub %in% Outliers_mm2[Outliers_mm2$category == "SC",]$Sub & MM2_Output$category == "SC", ]$mm2_z = NA #2 first removed (two below zero)
BiometricThin_Pheno = merge(BioMetricThin,MM2_Output, by = c("Sub","category"), all = TRUE)
#prepare data for the univariate Biometric model fitting
BiometricThinWide_Pheno = BiometricThin_Pheno%>%
  select(-Sub)%>%
  pivot_wider(names_from = "SibId", values_from = c("avgP","varP","mm2_r","mm2_z"))%>%
  mutate(Zygosity = ifelse(Zygosity == "MZ", 1,3))
#save
write_csv(BiometricThinWide_Pheno,sprintf("%s/%s/01_Germine_2015/05_twin_univariate.csv",wdOA,wdOA_output))

#MULTIVARIATE CTD####
#residualize age, sex
BiometricThin_Pheno_avgP = BiometricThin_Pheno%>%
  select(Sub, FamId, Age, Sex, avgP,Zygosity,category, SibId)%>%
  pivot_wider(names_from = category, values_from = c("avgP"))
#residualize age, sex
BiometricThin_Pheno_avgP_res = umx_residualize(c("AO","FA_TOT","SC"), c("Sex", "Age"),data = BiometricThin_Pheno_avgP)
#normalize for multivariate modeling
BiometricThin_Pheno_avgP_res = BiometricThin_Pheno_avgP_res %>% as.data.frame() %>%  mutate(AO = ((AO-mean(AO, na.rm = T))/sd(AO,na.rm = T)), 
                                                                                            FA_TOT = ((FA_TOT-mean(FA_TOT,na.rm = T))/sd(FA_TOT,na.rm = T)),
                                                                                            SC = ((SC-mean(SC,na.rm = T))/sd(SC,na.rm = T))
                                                                                            )
#residualize age, sex
BiometricThin_Pheno_mm2z = BiometricThin_Pheno%>%
  select(Sub,FamId, Age, Sex, mm2_z,Zygosity,category, SibId)%>%
  pivot_wider(names_from = category, values_from = c("mm2_z"))
#residualize age, sex
BiometricThin_Pheno_mm2z_res = umx_residualize(c("AO","FA_TOT","SC"), c("Sex", "Age"),data = BiometricThin_Pheno_mm2z)
#normalize for multivariate modeling
BiometricThin_Pheno_mm2z_res = BiometricThin_Pheno_mm2z_res %>% as.data.frame() %>%  mutate(AO = ((AO-mean(AO, na.rm = T))/sd(AO,na.rm = T)), 
                                                                                            FA_TOT = ((FA_TOT-mean(FA_TOT,na.rm = T))/sd(FA_TOT,na.rm = T)),
                                                                                            SC = ((SC-mean(SC,na.rm = T))/sd(SC,na.rm = T))
)
#prepare data for the multivariate Biometric model fitting
BiometricThinWide_Pheno_avgP_res = BiometricThin_Pheno_avgP_res%>%
  select(-Sub)%>%
  pivot_wider(names_from = "SibId", values_from = c("AO","FA_TOT","SC"))%>%
  mutate(Zygosity = ifelse(Zygosity == "MZ", 1,3))

BiometricThinWide_Pheno_mm2_res = BiometricThin_Pheno_mm2z_res%>%
  select(-Sub)%>%
  pivot_wider(names_from = "SibId", values_from = c("AO","FA_TOT","SC"))%>%
  mutate(Zygosity = ifelse(Zygosity == "MZ", 1,3))

#save
write_csv(BiometricThinWide_Pheno_avgP_res,sprintf("%s/%s/01_Germine_2015/05_twin_eb_multivariate.csv",wdOA,wdOA_output))
write_csv(BiometricThinWide_Pheno_mm2_res,sprintf("%s/%s/01_Germine_2015/05_twin_tt_multivariate.csv",wdOA,wdOA_output))

