#Author: Giacomo Bignardi
#Adapted from: NA
#Last modified: 2023-10-10
#
#
#
#Description:
#Program: twinDfs ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)
library(umx)
library(patchwork)


#clean working environment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

#load dataFrames:
BioMetric  = read_csv(sprintf("%s/%s/01_Germine_2015/01_CTD_Germine2015.csv", wdOA,wdOA_output))
BioMetric_val  = read_csv(sprintf("%s/%s/02_Sutherland_2020/01_CTD_Sutherland2020.csv", wdOA,wdOA_output))

#DISCOVERY####
#apply exclusion criteria
#average repeated measure first
Fac_BioMetric_avg1 = BioMetric%>%
  filter(category == "FA_TOT", intraR_FA >=.5)%>%
  filter(is_repeated == 1,SibId == 1)%>% 
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Fac_BioMetric1 = Fac_BioMetric_avg1 %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Sce_BioMetric_avg1  = BioMetric%>%
  filter(category == "SC", intraR_SC >=.5)%>%
  filter(is_repeated == 1,SibId == 1)%>% 
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Sce_BioMetric1 = Sce_BioMetric_avg1 %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Obj_BioMetric_avg1 = BioMetric%>%
  filter(category == "AO", intraR_AO >=.5)%>%
  filter(is_repeated == 1,SibId == 1)%>% 
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Obj_BioMetric1 = Obj_BioMetric_avg1 %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Fac_BioMetric1 = Fac_BioMetric1%>%
  mutate(category = "faces")
Sce_BioMetric1 = Sce_BioMetric1%>%
  mutate(category = "scenes")
Obj_BioMetric1 = Obj_BioMetric1%>%
  mutate(category = "abstracts")

#repeated image
#apply exclusion criteria
#average repeated measure first
Fac_BioMetric_avg2 = BioMetric%>%
  filter(category == "FA_TOT", intraR_FA >=.5)%>%
  filter(is_repeated == 2,SibId == 1)%>% 
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Fac_BioMetric2 = Fac_BioMetric_avg2 %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Sce_BioMetric_avg2  = BioMetric%>%
  filter(category == "SC", intraR_SC >=.5)%>%
  filter(is_repeated == 2,SibId == 1)%>% 
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Sce_BioMetric2 = Sce_BioMetric_avg2 %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Obj_BioMetric_avg2 = BioMetric%>%
  filter(category == "AO", intraR_AO >=.5)%>%
  filter(is_repeated == 2,SibId == 1)%>% 
  group_by(FamId_2,Item_2,FamId,SibId,Sex,Age,Zygosity,category,domain_2)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Obj_BioMetric2 = Obj_BioMetric_avg2 %>% 
  group_by(FamId_2,FamId,SibId,Sex,Age,Zygosity,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "FamId_2")

Fac_BioMetric2 = Fac_BioMetric2%>%
  mutate(category = "faces")
Sce_BioMetric2 = Sce_BioMetric2%>%
  mutate(category = "scenes")
Obj_BioMetric2 = Obj_BioMetric2%>%
  mutate(category = "abstracts")

Fac_BioMetric = merge(Fac_BioMetric1,Fac_BioMetric2, by = c("Sub","FamId","SibId", "Sex","Age", "Zygosity","category"))
Sce_BioMetric = merge(Sce_BioMetric1,Sce_BioMetric2, by = c("Sub","FamId","SibId", "Sex","Age", "Zygosity","category"))
Abs_BioMetric = merge(Obj_BioMetric1,Obj_BioMetric2, by = c("Sub","FamId","SibId", "Sex","Age", "Zygosity","category"))

#calculate test-retest
r_eb_faces_Germine = cor.test(as.numeric(Fac_BioMetric$avgP.x),as.numeric(Fac_BioMetric$avgP.y))
r_eb_scenes_Germine = cor.test(as.numeric(Sce_BioMetric$avgP.x),as.numeric(Sce_BioMetric$avgP.y))
r_eb_abstracts_Germine = cor.test(as.numeric(Abs_BioMetric$avgP.x),as.numeric(Abs_BioMetric$avgP.y))

#calculate test-retest icc 2,k (reliability using the average of the two)
icc_eb_faces_Germine = psych::ICC(Fac_BioMetric %>% select(avgP.x,avgP.y))
icc_eb_scenes_Germine = psych::ICC(Sce_BioMetric %>% select(avgP.x,avgP.y))
icc_eb_abstracts_Germine = psych::ICC(Abs_BioMetric %>% select(avgP.x,avgP.y))

#VALIDATION####
#apply exclusion criteria
#average repeated measure first
Fac_BioMetric_avg1_val = BioMetric_val%>%
  filter(category == "FA", intraR_FA >=.5)%>%
  filter(is_repeated == 1,twinN == "A")%>% 
  select(-ZygosityStatus)%>%
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Fac_BioMetric1_val = Fac_BioMetric_avg1_val %>% 
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")

Sce_BioMetric_avg1_val  = BioMetric_val%>%
  filter(category == "SC", intraR_SC >=.5)%>%
  filter(is_repeated == 1,twinN == "A")%>% 
  select(-ZygosityStatus)%>%
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Sce_BioMetric1_val = Sce_BioMetric_avg1_val %>% 
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")



Fac_BioMetric1_val = Fac_BioMetric1_val%>%
  mutate(category = "faces")
Sce_BioMetric1_val = Sce_BioMetric1_val%>%
  mutate(category = "scenes")

#repeated image
#apply exclusion criteria
#average repeated measure first
Fac_BioMetric_avg2_val = BioMetric_val%>%
  filter(category == "FA", intraR_FA >=.5)%>%
  filter(is_repeated ==2,twinN == "A")%>% 
  select(-ZygosityStatus)%>%
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Fac_BioMetric2_val = Fac_BioMetric_avg2_val %>% 
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")

Sce_BioMetric_avg2_val  = BioMetric_val%>%
  filter(category == "SC", intraR_SC >=.5)%>%
  filter(is_repeated == 2,twinN == "A")%>% 
  select(-ZygosityStatus)%>%
  group_by(refCode,Item,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(Value = mean(Value))%>%
  ungroup()
Sce_BioMetric2_val = Sce_BioMetric_avg2_val %>% 
  group_by(refCode,twinN,Sex,AgeAtTestingCalculated,TwinType,category)%>%
  summarise(avgP =mean(Value), varP = var(Value))%>%
  ungroup()%>%
  rename(Sub = "refCode")

Fac_BioMetric2_val = Fac_BioMetric2_val%>%
  mutate(category = "faces")
Sce_BioMetric2_val = Sce_BioMetric2_val%>%
  mutate(category = "scenes")

Fac_BioMetric_val = merge(Fac_BioMetric1_val,Fac_BioMetric2_val, by = c("Sub","twinN","Sex", "AgeAtTestingCalculated","TwinType", "category"))
Sce_BioMetric_val = merge(Sce_BioMetric1_val,Sce_BioMetric2_val, by = c("Sub","twinN","Sex", "AgeAtTestingCalculated","TwinType", "category"))

#calculate test-retest
r_eb_faces_Sutherland = cor.test(as.numeric(Fac_BioMetric_val$avgP.x),as.numeric(Fac_BioMetric_val$avgP.y))
r_eb_scenes_Sutherland = cor.test(as.numeric(Sce_BioMetric_val$avgP.x),as.numeric(Sce_BioMetric_val$avgP.y))

#calculate test-retest icc 2,k (reliability using the average of the two)
icc_eb_faces_Sutherland = psych::ICC(Fac_BioMetric_val %>% select(avgP.x,avgP.y))
icc_eb_scenes_Sutherland = psych::ICC(Sce_BioMetric_val %>% select(avgP.x,avgP.y))

#SAVE####
save(r_eb_faces_Germine, r_eb_scenes_Germine, r_eb_abstracts_Germine,
       icc_eb_faces_Germine,icc_eb_scenes_Germine, icc_eb_abstracts_Germine,
       r_eb_faces_Sutherland,r_eb_scenes_Sutherland,
       icc_eb_faces_Sutherland,icc_eb_scenes_Sutherland, file = sprintf("%s/%s/04_review/02_report_withinday_reliability_eb.RData",wdOA,wdOA_output))#save output
