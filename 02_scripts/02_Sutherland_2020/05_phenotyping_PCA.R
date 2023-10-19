#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description: repeat procedure in highlighted in 01_Germine_2015/05_phenotyping_PCA.R
# additional step: residualize for control bias and typicality
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
BioMetric  = read_csv(sprintf("%s/%s/02_Sutherland_2020/01_CTD_Sutherland2020.csv", wdOA,wdOA_output))
Fac_MM2_Output = read_csv(sprintf("%s/%s/02_Sutherland_2020/04_taste_typicality_faces_Sutherland2020.csv", wdOA,wdOA_output))
Sce_MM2_Output = read_csv(sprintf("%s/%s/02_Sutherland_2020/04_taste_typicality_scenes_Sutherland2020.csv", wdOA,wdOA_output))
Con_MM2_Output = read_csv(sprintf("%s/%s/02_Sutherland_2020/04_taste_typicality_control_Sutherland2020.csv", wdOA,wdOA_output))

#EVALUATION-BIAS#####
#Compute average per participant after standardization
####_faces####
#prepare dataframe (apply exclusion criteria)
Fac_BioMetricThin= BioMetric%>%
  select(-ZygosityStatus)%>%
  filter(category == "FA", intraR_FA >=.5)%>% #exclusion criteria
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))# average across repeated measure first
#calculate overall pleasantess
Fac_BioMetricSlim_avg = Fac_BioMetricThin%>%
  ungroup()%>%
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")

####_scenes####
Sce_BioMetricThin = BioMetric%>%
  select(-ZygosityStatus)%>%
  filter(category == "SC", intraR_SC >=.5)%>%
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))
#calculate overall pleasantess
Sce_BioMetricSlim_avg = Sce_BioMetricThin%>%
  ungroup()%>%
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")

####_Control####
Con_BioMetricThin = BioMetric%>%
  select(-ZygosityStatus)%>%
  filter(category == "DO")%>%
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))
#calculate overall dominance
Con_BioMetricSlim_avg = Con_BioMetricThin%>%
  ungroup()%>%
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")

#bind togheter 
BioMetricThin_avg = rbind(Fac_BioMetricSlim_avg,Sce_BioMetricSlim_avg,Con_BioMetricSlim_avg);nrow(BioMetricThin_avg)#check they are the same

####MM2####
Fac_MM2_Output = Fac_MM2_Output%>%
  mutate(category = "FA")
Sce_MM2_Output = Sce_MM2_Output%>%
  mutate(category = "SC")
Con_MM2_Output = Con_MM2_Output%>%
  mutate(category = "DO")
MM2_Output = rbind(Sce_MM2_Output,Fac_MM2_Output,Con_MM2_Output)


#####Outliers####
## Residualization and PCA are both sensitive to extreme outliers: as such remove extreme outliers for specific domain from subsequent analysis
#exclude extreme outliers
Outliers_mm2 = MM2_Output%>%
  mutate(twinN = substr(Sub,nchar(Sub),nchar(Sub)))   %>% 
  group_by(category, twinN)%>%
  rstatix::identify_outliers(mm2_z)%>%
  filter(is.extreme == T) %>% 
 select(category,Sub)
Outliers_mm2_FA = Outliers_mm2 %>% filter(category == "FA") %>% pull(Sub)#4

Outliers_p = BioMetricThin_avg%>%
  group_by(category, twinN)%>%
  rstatix::identify_outliers(avgP)%>%
  filter(is.extreme == T) %>% 
  filter(category != "DO") %>% 
  select(category,Sub)
Outliers_p_SC = Outliers_p %>% filter(category == "SC") %>% pull(Sub)#2

#PCA####
set.seed(123)
#####compute PC####
Fac_BioMetric_pc_t1 = Fac_BioMetricThin %>% 
  ungroup() %>% 
  filter(twinN == "A") %>% 
  filter(!( refCode %in% c(Outliers_mm2_FA))) %>% 
  select( refCode, Value ,Item) %>% 
  pivot_wider(values_from ="Value", names_from = "Item") #prepare wide dataframe (row individuals, column images)
Fac_PC_t1 = prcomp(Fac_BioMetric_pc_t1[,-1],robust = c("S-estimator")) #Principal Componet Analysis
Fac_PCscore_t1 = cbind(Fac_BioMetric_pc_t1[,1],Fac_PC_t1$x) #Extract score per individual and bind with individual ID
Sce_BioMetric_pc_t1 = Sce_BioMetricThin %>% 
  ungroup() %>% 
  filter(twinN == "A") %>% 
  filter(!( refCode %in% c(Outliers_p_SC))) %>% 
  select( refCode, Value ,Item) %>% 
  pivot_wider(values_from ="Value", names_from = "Item")
Sce_PC_t1 = prcomp(Sce_BioMetric_pc_t1[,-1])
Sce_PCscore_t1 = cbind(Sce_BioMetric_pc_t1[,1],Sce_PC_t1$x)

Fac_BioMetric_pc_t2 = Fac_BioMetricThin %>% 
  ungroup() %>% 
  filter(twinN == "B") %>% 
  filter(!( refCode %in% c(Outliers_mm2_FA))) %>% 
  select( refCode, Value ,Item) %>% 
  pivot_wider(values_from ="Value", names_from = "Item") #prepare wide dataframe (row individuals, column images)
Fac_PC_t2 = prcomp(Fac_BioMetric_pc_t2[,-1]) #Principal Componet Analysis
Fac_PCscore_t2 = cbind(Fac_BioMetric_pc_t2[,1],Fac_PC_t2$x) #Extract score per individual and bind with individual ID
Sce_BioMetric_pc_t2 = Sce_BioMetricThin %>% 
  ungroup() %>% 
  filter(twinN == "B") %>% 
  filter(!( refCode %in% c(Outliers_p_SC))) %>% 
  select( refCode, Value ,Item) %>% 
  pivot_wider(values_from ="Value", names_from = "Item")
Sce_PC_t2 = prcomp(Sce_BioMetric_pc_t2[,-1])
Sce_PCscore_t2 = cbind(Sce_BioMetric_pc_t2[,1],Sce_PC_t2$x)

#####Proportion of variance####
#Variance explained
FA_PC_t1  = Fac_PC_t1$sdev^2 / sum(Fac_PC_t1$sdev^2)
SC_PC_t1  = Sce_PC_t1$sdev^2 / sum(Sce_PC_t1$sdev^2)

FA_PC_t2  = Fac_PC_t2$sdev^2 / sum(Fac_PC_t2$sdev^2)
SC_PC_t2  = Sce_PC_t2$sdev^2 / sum(Sce_PC_t2$sdev^2)

#total variance
PC12_FA_var_explained_t1_Sutherland = sum(FA_PC_t1[1:2])
PC12_SC_var_explained_t1_Sutherland = sum(SC_PC_t1[1:2])

PC12_FA_var_explained_t2_Sutherland = sum(FA_PC_t2[1:2])
PC12_SC_var_explained_t2_Sutherland = sum(SC_PC_t2[1:2])

#####Correlation with facets####
#bind twin one and twin two and prepare to add info to dataframe
Fac_PCscore = rbind(Fac_PCscore_t1[,1:3], Fac_PCscore_t2[,1:3])
Sce_PCscore = rbind(Sce_PCscore_t1[,1:3], Sce_PCscore_t2[,1:3])

Fac_BioMetric_PC_avg = merge(Fac_BioMetricSlim_avg %>%  filter(!( Sub %in% Outliers_mm2_FA)) ,Fac_PCscore %>% rename(Sub   = "refCode"), by = c("Sub"), all = T)
Sce_BioMetric_PC_avg = merge(Sce_BioMetricSlim_avg %>%  filter(!( Sub %in% Outliers_p_SC)),Sce_PCscore %>% rename(Sub   = "refCode"), by = c("Sub"), all = T)

Fac_BioMetric_PC_mm2 = merge(Fac_MM2_Output %>% select(-category) %>%  filter(!( Sub %in% Outliers_mm2_FA)),Fac_PCscore %>% rename(Sub   = "refCode"), by = c("Sub"), all = T);Fac_BioMetric_PC_mm2 =  Fac_BioMetric_PC_mm2 %>% mutate(twinN = substr(Sub,nchar(Sub),nchar(Sub)))
Sce_BioMetric_PC_mm2 = merge(Sce_MM2_Output %>% select(-category) %>% filter(!( Sub %in% Outliers_p_SC)),Sce_PCscore %>% rename(Sub   = "refCode"), by = c("Sub"), all = T);Sce_BioMetric_PC_mm2 =  Sce_BioMetric_PC_mm2 %>% mutate(twinN = substr(Sub,nchar(Sub),nchar(Sub)))


#Correlation between PCs and Facets of aesthetic value
PC_FA_cor_eb_t1_Sutherland = cor.test(-(Fac_BioMetric_PC_avg %>% filter(twinN == "A") %>% pull(PC1)),  Fac_BioMetric_PC_avg %>% filter(twinN == "A") %>% pull(avgP))#flipped
PC_SC_cor_eb_t1_Sutherland = cor.test(Sce_BioMetric_PC_avg %>% filter(twinN == "A") %>% pull(PC1),  Sce_BioMetric_PC_avg %>% filter(twinN == "A") %>% pull(avgP))
PC_FA_cor_eb_t2_Sutherland = cor.test(Fac_BioMetric_PC_avg %>% filter(twinN == "B") %>% pull(PC1),  Fac_BioMetric_PC_avg %>% filter(twinN == "B") %>% pull(avgP))
PC_SC_cor_eb_t2_Sutherland = cor.test(Sce_BioMetric_PC_avg %>% filter(twinN == "B") %>% pull(PC1),  Sce_BioMetric_PC_avg %>% filter(twinN == "B") %>% pull(avgP))

PC_FA_cor_tt_t1_Sutherland = cor.test(-(Fac_BioMetric_PC_mm2 %>% filter(twinN == "A") %>% pull(PC2)),  Fac_BioMetric_PC_mm2 %>% filter(twinN == "A") %>% pull(mm2_z))#flipped
PC_SC_cor_tt_t1_Sutherland = cor.test(-(Sce_BioMetric_PC_mm2 %>% filter(twinN == "A") %>% pull(PC2)),  Sce_BioMetric_PC_mm2 %>% filter(twinN == "A") %>% pull(mm2_z))#flipped
PC_FA_cor_tt_t2_Sutherland = cor.test(Fac_BioMetric_PC_mm2 %>% filter(twinN == "B") %>% pull(PC2),  Fac_BioMetric_PC_mm2 %>% filter(twinN == "B") %>% pull(mm2_z))
PC_SC_cor_tt_t2_Sutherland = cor.test(Sce_BioMetric_PC_mm2 %>% filter(twinN == "B") %>% pull(PC2),  Sce_BioMetric_PC_mm2 %>% filter(twinN == "B") %>% pull(mm2_z))

##REPORT:PCA####
save(PC_FA_cor_eb_t1_Sutherland, PC_SC_cor_eb_t1_Sutherland,
     PC_FA_cor_eb_t2_Sutherland, PC_SC_cor_eb_t2_Sutherland,
     PC_FA_cor_tt_t1_Sutherland, PC_SC_cor_tt_t1_Sutherland,
     PC_FA_cor_tt_t2_Sutherland, PC_SC_cor_tt_t2_Sutherland,
     PC12_FA_var_explained_t1_Sutherland,
     PC12_SC_var_explained_t1_Sutherland,
     PC12_FA_var_explained_t2_Sutherland,
     PC12_SC_var_explained_t2_Sutherland,
     file = sprintf("%s/%s/02_Sutherland_2020/05_report_PCA_val.Rdata",wdOA,wdOA_output))


#####F1e####
p_Fac_BioMetric_PC_avg_t1_Sut2020 = Fac_BioMetric_PC_avg %>% 
  filter(twinN == "A") %>% 
  ggplot(aes(y = -PC1,avgP)) + 
  xlim(1,9)+
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
    y = "(-)",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Sce_BioMetric_PC_avg_t1_Sut2020 = Sce_BioMetric_PC_avg %>% 
  filter(twinN == "A") %>% 
  ggplot(aes(y =PC1,avgP)) + 
  xlim(1,9)+
  geom_point(fill = "black",colour = "white",pch=21,  size = 1.5) + 
  labs(#title = "Scenes",
    y = "",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Fac_BioMetric_PC_mm2_t1_Sut2020 = Fac_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(twinN == "A") %>% 
  ggplot(aes(y =-PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
    y = "(-)",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Sce_BioMetric_PC_mm2_t1_Sut2020 = Sce_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(twinN == "A") %>% 
  ggplot(aes(y =-PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Scenes",
    y = "(-)",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Fac_BioMetric_PC_avg_t2_Sut2020 = Fac_BioMetric_PC_avg %>% 
  filter(twinN == "B") %>% 
  ggplot(aes(y = PC1,avgP)) + 
  xlim(1,9)+
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
    y = "",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Sce_BioMetric_PC_avg_t2_Sut2020 = Sce_BioMetric_PC_avg %>% 
  filter(twinN == "B") %>% 
  ggplot(aes(y =PC1,avgP)) + 
  xlim(1,9)+
  geom_point(fill = "black",colour = "white",pch=21,  size = 1.5) + 
  labs(#title = "Scenes",
    y = "",
    x = "Evaluation-bias") +
  theme_classic(base_size = 16)

p_Fac_BioMetric_PC_mm2_t2_Sut2020 = Fac_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(twinN == "B") %>% 
  ggplot(aes(y =PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Faces",
    y = "",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 

p_Sce_BioMetric_PC_mm2_t2_Sut2020 = Sce_BioMetric_PC_mm2 %>% 
  mutate(SibId = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(twinN == "B") %>% 
  ggplot(aes(y =PC2,mm2_z)) + 
  geom_point(fill = "black",colour = "white",pch=21, size = 1.5) + 
  labs(#title = "Scenes",
    y = "",
    x = "taste-typicality") +
  theme_classic(base_size = 16) 


#Save pc
save(p_Fac_BioMetric_PC_avg_t1_Sut2020, p_Fac_BioMetric_PC_avg_t2_Sut2020, p_Sce_BioMetric_PC_avg_t1_Sut2020, p_Sce_BioMetric_PC_avg_t2_Sut2020,
     p_Fac_BioMetric_PC_mm2_t1_Sut2020, p_Fac_BioMetric_PC_mm2_t2_Sut2020, p_Sce_BioMetric_PC_mm2_t1_Sut2020, p_Sce_BioMetric_PC_mm2_t2_Sut2020, file = 
     sprintf("%s/%s/02_Sutherland_2020/05_sup_S5c_PC_Sutherland2020.Rdata",wdOA,wdOA_output))



#####FS4####
#Plot PC1 PC2 score for abstract and overlay to pleasantess and taste-typicality scores
p_PC1_FA_sut2020 = 
  ggplot(Fac_BioMetric_PC_avg %>% filter(twinN == "A"), aes(-PC1, -PC2, color = avgP)) +
  geom_point() +
  labs(x = paste0("(-) PC1 (",round(FA_PC_t1[1],2),"%)"),
       y = paste0("(-) PC2 (",round(FA_PC_t1[2],2),"%)"),
       color = "Evaluation-bias",
       subtitle = "Faces")+
  scale_color_viridis_c(option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC2_FA_sut2020 =ggplot(Fac_BioMetric_PC_mm2 %>% filter(twinN == "A"), aes(-PC1, -PC2, color = mm2_z)) +
  geom_point() +
  labs(x = paste0("(-) PC1 (",round(FA_PC_t1[1],2),"%)"),
       y =  paste0(" (-)PC2 (",round(FA_PC_t1[2],2),"%)"),
       color = "Taste-typicality",
       subtitle = "Faces")+
  scale_color_viridis_c(option = "cividis") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC1_SC_sut2020 = 
  ggplot(Sce_BioMetric_PC_avg %>% filter(twinN == "A"), aes(PC1, -PC2, color = avgP)) +
  geom_point() +
  labs(x = paste0("PC1 (",round(SC_PC_t1[1],2),"%)"),
       y = paste0("(-) PC2 (",round(SC_PC_t1[2],2),"%)"),
       color = "Evaluation-bias",
       subtitle = "Scenes")+
  scale_color_viridis_c(option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

p_PC2_SC_sut2020 =ggplot(Sce_BioMetric_PC_mm2 %>% filter(twinN == "A"), aes(PC1, -PC2, color = mm2_z)) +
  geom_point() +
  labs(x = paste0("PC1 (",round(SC_PC_t1[1],2),"%)"),
       y =  paste0("(-) PC2 (",round(SC_PC_t1[2],2),"%)"),
       color = "Taste-typicality",
       subtitle = "Scenes")+
  scale_color_viridis_c(option = "cividis") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  theme_minimal(base_size = 18)

#scree plot
#take variance explained from each component (select first 10)
PC_t1_sp = rbind(
  FA_PC_t1_sp = data.frame(prop_variance = FA_PC_t1[1:10], component = rep(1:10), visual_domain = "faces"),
  SC_PC_t1_sp = data.frame(prop_variance = SC_PC_t1[1:10],component =  rep(1:10),  visual_domain = "scenes")
)

PC_t1_sp = PC_t1_sp %>% mutate(visual_domain = factor(visual_domain, levels= c("faces","scenes")))

#take variance explained from each component (select first 10) for twin 2
PC_t2_sp = rbind(
  FA_PC_t2_sp = data.frame(prop_variance = FA_PC_t2[1:10], component = rep(1:10), visual_domain = "faces"),
  SC_PC_t2_sp = data.frame(prop_variance = SC_PC_t2[1:10],component =  rep(1:10),  visual_domain = "scenes")
)

PC_t2_sp = PC_t2_sp %>% mutate(visual_domain = factor(visual_domain, levels= c("faces","scenes")))

#Scree plot
p_sp_1_sut2020 = ggplot(data=PC_t1_sp, aes(x=component, y=prop_variance, color = visual_domain)) +
  geom_line()+
  scale_x_continuous(breaks=seq(1, 10,1))+
  ylim(c(0,0.5))+
  scale_color_brewer(palette = "Dark2")+
  geom_point() +
  theme_classic()
p_sp_2_sut2020 = ggplot(data=PC_t2_sp, aes(x=component, y=prop_variance, color = visual_domain)) +
  geom_line()+
  scale_x_continuous(breaks=seq(1, 10,1))+
  ylim(c(0,0.5))+
  scale_color_brewer(palette = "Dark2")+
  geom_point() +
  theme_classic()

save(p_PC2_SC_sut2020,
     p_PC1_SC_sut2020,
     p_PC2_FA_sut2020, 
     p_PC1_FA_sut2020,
     p_sp_1_sut2020,
     p_sp_2_sut2020,
     file = sprintf("%s/%s/02_Sutherland_2020/05_sup_FS5_screeplot_PCA_sutherland2020.Rdata",wdOA,wdOA_output))


#OUTLIERS####
#convert to NA for later analysis (NOTE: to obtain results without outlier removal just uncomment this section)
BioMetricThin_avg[BioMetricThin_avg$Sub %in% Outliers_p[Outliers_p$category == "FA",]$Sub & BioMetricThin_avg$category == "FA", ]$avgP = NA #0 removed
BioMetricThin_avg[BioMetricThin_avg$Sub %in% Outliers_p[Outliers_p$category == "SC",]$Sub & BioMetricThin_avg$category == "SC", ]$avgP = NA #1 removed

#convert to NA (mantain pair member in biometric model)
MM2_Output[MM2_Output$Sub %in% Outliers_mm2[Outliers_mm2$category == "FA",]$Sub & MM2_Output$category == "FA", ]$mm2_z = NA #4 removed
MM2_Output[MM2_Output$Sub %in% Outliers_mm2[Outliers_mm2$category == "SC",]$Sub & MM2_Output$category == "SC", ]$mm2_z = NA #0 removed

#prepare chr as required by SAT model fitting
BioMetricThin_avg = BioMetricThin_avg%>%mutate(Sub = substr(Sub,1,nchar(Sub)-2))
MM2_Output = MM2_Output%>%mutate(twinN = substr(Sub,nchar(Sub),nchar(Sub)),
                                 Sub = substr(Sub,1,nchar(Sub)-2))

BiometricThin_Pheno = merge(BioMetricThin_avg,MM2_Output, by = c("Sub","twinN","category"), all = TRUE)

####RESIDUALIZATION####
#regress bias for P
BiometricThin_Pheno_avgP = BiometricThin_Pheno%>%
  select(Sub, AgeAtTestingCalculated, Sex, avgP,TwinType,category, twinN)%>%
  pivot_wider(names_from = category, values_from = c("avgP"))
#residualize per bias
BiometricThin_Pheno_avgP_res = umx_residualize(c("FA","SC"), c("DO"),data = BiometricThin_Pheno_avgP) #dominance
BiometricThin_Pheno_avgP_res_SexAge = umx_residualize(c("FA","SC"), c("AgeAtTestingCalculated", "Sex"),data = BiometricThin_Pheno_avgP) #sex and age
BiometricThin_Pheno_avgP_res_multive = umx_residualize(c("FA","SC"), c("AgeAtTestingCalculated", "Sex"),data = BiometricThin_Pheno_avgP_res) #both
#normalize for multivariate modeling
BiometricThin_Pheno_avgP_res_SexAge = BiometricThin_Pheno_avgP_res_SexAge %>% as.data.frame() %>%  mutate(FA = ((FA-mean(FA, na.rm = T))/sd(FA,na.rm = T)), SC = ((SC-mean(SC,na.rm = T))/sd(SC,na.rm = T)))
BiometricThin_Pheno_avgP_res_multive = BiometricThin_Pheno_avgP_res_multive %>% as.data.frame() %>%   mutate(FA = ((FA-mean(FA, na.rm = T))/sd(FA,na.rm = T)), SC = ((SC-mean(SC,na.rm = T))/sd(SC,na.rm = T)))


p_avgP_faces_val = ggplot(BiometricThin_Pheno_avgP) + 
  geom_histogram(bins = 50, aes(x = FA,), position="identity", colour='white') + 
  labs(x = "pleasanteness faces val")+
  scale_y_continuous(limits= c(0,90),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(-1,0,1))+
  theme_classic(base_size=16)

p_avgP_scenes_val =ggplot(BiometricThin_Pheno_avgP) + 
  geom_histogram(bins = 50, aes(x = SC), position="identity", colour='white') + 
  labs(x = "pleasanteness")+
  scale_y_continuous(limits= c(0,90),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(-1,0,1))+
  theme_classic(base_size=16)

p_avgP_faces_val
p_avgP_scenes_val

#Pregress bias for mm2
BiometricThin_Pheno_mm2 = BiometricThin_Pheno%>%
  select(Sub, AgeAtTestingCalculated, Sex, mm2_z,TwinType,category, twinN)%>%
  pivot_wider(names_from = category, values_from = c("mm2_z"))
#residualize per bias
BiometricThin_Pheno_mm2_res = umx_residualize(c("FA","SC"), c("DO"),data = BiometricThin_Pheno_mm2)
BiometricThin_Pheno_mm2_res_SexAge = umx_residualize(c("FA","SC"), c("AgeAtTestingCalculated", "Sex"),data = BiometricThin_Pheno_mm2)
BiometricThin_Pheno_mm2_res_multive = umx_residualize(c("FA","SC"), c("AgeAtTestingCalculated", "Sex"),data = BiometricThin_Pheno_mm2_res)

scale(BiometricThin_Pheno_mm2_res_SexAge$FA)
#normalize for multivariate modeling
BiometricThin_Pheno_mm2_res_SexAge = BiometricThin_Pheno_mm2_res_SexAge %>% as.data.frame() %>%  mutate(FA = ((FA-mean(FA, na.rm = T))/sd(FA,na.rm = T)), SC = ((SC-mean(SC,na.rm = T))/sd(SC,na.rm = T)))
BiometricThin_Pheno_mm2_res_multive = BiometricThin_Pheno_mm2_res_multive %>% as.data.frame() %>% mutate(FA = ((FA-mean(FA, na.rm = T))/sd(FA,na.rm = T)), SC = ((SC-mean(SC,na.rm = T))/sd(SC,na.rm = T)))

p_mm2_faces_val = ggplot(BiometricThin_Pheno_mm2) + 
  geom_histogram(bins = 50, aes(x = FA), position="identity", colour='white') + 
  labs(x = "t-t faces val.")+
  scale_y_continuous(limits= c(0,90),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(-1,0,1))+
  theme_classic(base_size=16)

Outliers_mm2[Outliers_mm2$category == "FA",]$Sub
p_mm2_scenes_val = ggplot(BiometricThin_Pheno_mm2) + 
  geom_histogram(bins = 50, aes(x = SC), position="identity", colour='white') + 
  labs(x = "t-t scenes val.")+
  scale_y_continuous(limits= c(0,90),breaks = c(0,40,80))+
  scale_x_continuous(breaks = c(-1,0,1))+
  theme_classic(base_size=16)

#save for later plotting
p_mm2_faces_val
p_mm2_scenes_val

#save for later phenotypic correlations
save(BiometricThin_Pheno_avgP_res_SexAge,BiometricThin_Pheno_avgP_res_multive,BiometricThin_Pheno_mm2_res_SexAge,BiometricThin_Pheno_mm2_res_multive,file = sprintf("%s/%s/02_Sutherland_2020/05_phenotypic_cor.Rdata",wdOA,wdOA_output))

#reshape to wide for univariate  and multivariate modeling
BiometricThinWide_Pheno_avgP= BiometricThin_Pheno_avgP%>% #no residuals
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))
BiometricThinWide_Pheno_avgP_res = BiometricThin_Pheno_avgP_res%>% #residuals (only dominance)
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))
BiometricThinWide_Pheno_avgP_res_SexAge = BiometricThin_Pheno_avgP_res_SexAge%>% #residuals (only sex and age)
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))
BiometricThinWide_Pheno_avgP_res_multive = BiometricThin_Pheno_avgP_res_multive%>% #residuals (dominance, sex and age)
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))

BiometricThinWide_Pheno_mm2 = BiometricThin_Pheno_mm2%>%
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))
BiometricThinWide_Pheno_mm2_res = BiometricThin_Pheno_mm2_res%>%
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))
BiometricThinWide_Pheno_mm2_res_SexAge = BiometricThin_Pheno_mm2_res_SexAge%>%
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))
BiometricThinWide_Pheno_mm2_res_multive = BiometricThin_Pheno_mm2_res_multive%>%
  select(-c(DO))%>%
  pivot_wider(names_from = "twinN", values_from = c("FA","SC"))%>%
  mutate(TwinType = ifelse(TwinType == "MZ", 1,3))

#Vizualize residuals
ex_res = lm(SC~DO, data = BiometricThin_Pheno_avgP)
ex_res_agu = broom::augment(ex_res)
p_example_resid = ggplot(ex_res_agu, aes( DO, SC)) +
  labs(x = "Control-bias", y = "Evaluation-bias", title = "Scenes")+
  geom_segment(aes(xend = DO, yend = .fitted), color = "black", size = 0.5, alpha = 0.2, linetype = "dashed") +
  geom_point(alpha = .5, size = .75) +
  stat_smooth(method = lm) +
  theme_classic()

save(p_example_resid, file = sprintf("%s/%s/02_Sutherland_2020/05_sup_FS7_residuals_sutherland2020.Rdata",wdOA,wdOA_output))

#save
#save for univariate model:
#save non residualized scores
write_csv(BiometricThinWide_Pheno_mm2_res_multive,sprintf("%s/%s/02_Sutherland_2020/05_twin_multivariate_val.csv",wdOA,wdOA_output))
write_csv(BiometricThinWide_Pheno_avgP,sprintf("%s/%s/02_Sutherland_2020/05_twin_eb.csv",wdOA,wdOA_output))
write_csv(BiometricThinWide_Pheno_mm2,sprintf("%s/%s/02_Sutherland_2020/05_twin_tt.csv",wdOA,wdOA_output))
#save residualized scores
write_csv(BiometricThinWide_Pheno_avgP_res,sprintf("%s/%s/02_Sutherland_2020/05_twin_ebres.csv",wdOA,wdOA_output))
write_csv(BiometricThinWide_Pheno_mm2_res,sprintf("%s/%s/02_Sutherland_2020/05_twin_ttres.csv",wdOA,wdOA_output))
#save for multivariate model:
#save residualized scores (only sex and age)
write_csv(BiometricThinWide_Pheno_avgP_res_SexAge,sprintf("%s/%s/02_Sutherland_2020/05_twin_eb_multivariate.csv",wdOA,wdOA_output))
write_csv(BiometricThinWide_Pheno_mm2_res_SexAge,sprintf("%s/%s/02_Sutherland_2020/05_twin_tt_multivariate.csv",wdOA,wdOA_output))
#save residualized scores
write_csv(BiometricThinWide_Pheno_avgP_res_multive,sprintf("%s/%s/02_Sutherland_2020/05_twin_ebres_multivariate.csv",wdOA,wdOA_output))
write_csv(BiometricThinWide_Pheno_mm2_res_multive,sprintf("%s/%s/02_Sutherland_2020/05_twin_ttres_multivariate.csv",wdOA,wdOA_output))
