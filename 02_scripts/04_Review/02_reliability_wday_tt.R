#Author: Giacomo Bignardi
#Adapted from: NA
#Last modified: 2023-10-10
#
#
#
#Description:
#Program: xxx ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)

#clean working environment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

#load dataFrames:
MM2_TwinLong  = read_csv(sprintf("%s/%s/01_Germine_2015/01_MM2_Germine2015.csv", wdOA,wdOA_output))
MM2_TwinLong_val  = read_csv(sprintf("%s/%s/02_Sutherland_2020/01_MM2_Sutherland2020.csv", wdOA,wdOA_output))

#load functions:
source(sprintf("%s/%s/functions/MM2.R", wdOA,wdOA_scripts))

#DISCOVERY#####
MM2_TwinLong = as.data.frame(MM2_TwinLong)

#Prepare MM2_data for MM2 function
names(MM2_TwinLong) = c("Obj","Rating","Sub","Block","domain") #rename first 4 columns

MM2_TwinWide = MM2_TwinLong%>%
  pivot_wider(names_from = "Block", values_from = "Rating")%>%
  filter(!is.na(`2`))

Fac_MM2_TwinWide = as.data.frame(MM2_TwinWide%>%filter(domain == "FA_TOT")%>%
                                   select(-domain))
Sce_MM2_TwinWide = as.data.frame(MM2_TwinWide%>%filter(domain == "SC")%>%
                                   select(-domain)) 
Obj_MM2_TwinWide = as.data.frame(MM2_TwinWide%>%filter(domain == "AO")%>%
                                   select(-domain)) 
#extract info from the df
Fac_N = length(unique(Fac_MM2_TwinWide[,2])) #N participants
Fac_N_items = length(unique(Fac_MM2_TwinWide[,1]))

#extract info from the df
Sce_N = length(unique(Sce_MM2_TwinWide[,2])) #N participants
Sce_N_items = length(unique(Sce_MM2_TwinWide[,1]))

#extract info from the df
Obj_N = length(unique(Obj_MM2_TwinWide[,2])) #N participants
Obj_N_items = length(unique(Obj_MM2_TwinWide[,1]))

#prepare final matricies for images rated the first time
Fac_MM2_TwinMatrix1  = Fac_MM2_TwinWide %>%
  select(-`2`) %>% 
  pivot_wider(names_from = "Sub", values_from = "1")
Sce_MM2_TwinMatrix1  = Sce_MM2_TwinWide %>%
  select(-`2`) %>% 
  pivot_wider(names_from = "Sub", values_from = "1")
Obj_MM2_TwinMatrix1  = Obj_MM2_TwinWide%>%
  select(-`2`) %>% 
  pivot_wider(names_from = "Sub", values_from = "1")

####MM2####
#caclulate MM2 (Mean Muinus 2, see MM2.R function)
MM2_faces1 = MM2(Fac_MM2_TwinMatrix1); MM2_faces1$MM2
MM2_scenes1 = MM2(Sce_MM2_TwinMatrix1); MM2_scenes1$MM2
MM2_objects1 = MM2(Obj_MM2_TwinMatrix1); MM2_objects1$MM2

#prepare final matricies for images rated the second time
Fac_MM2_TwinMatrix2  = Fac_MM2_TwinWide %>%
  select(-`1`) %>% 
  pivot_wider(names_from = "Sub", values_from = "2")
Sce_MM2_TwinMatrix2  = Sce_MM2_TwinWide %>%
  select(-`1`) %>% 
  pivot_wider(names_from = "Sub", values_from = "2")
Obj_MM2_TwinMatrix2  = Obj_MM2_TwinWide%>%
  select(-`1`) %>% 
  pivot_wider(names_from = "Sub", values_from = "2")

####MM2####
#caclulate MM2 (Mean Muinus 2, see MM2.R function)
MM2_faces2 = MM2(Fac_MM2_TwinMatrix2); MM2_faces2$MM2
MM2_scenes2 = MM2(Sce_MM2_TwinMatrix2); MM2_scenes2$MM2
MM2_objects2 = MM2(Obj_MM2_TwinMatrix2); MM2_objects2$MM2

#merge mm2 stratified by repeated image
MM2_faces = merge(MM2_faces1$summary,MM2_faces2$summary, by = "Sub")
MM2_scenes = merge(MM2_scenes1$summary,MM2_scenes2$summary, by = "Sub")
MM2_objects = merge(MM2_objects1$summary,MM2_objects2$summary, by = "Sub")

#compute test retest on twin one only
MM2_faces = MM2_faces %>% as.data.frame() %>% 
  mutate(SibID = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibID == 1) %>% 
  select(-SibID)
MM2_scenes = MM2_scenes %>% as.data.frame() %>% 
  mutate(SibID = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibID == 1) %>% 
  select(-SibID)
MM2_objects = MM2_objects %>% as.data.frame() %>% 
  mutate(SibID = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibID == 1) %>% 
  select(-SibID)
#calculate test-retest
r_tt_faces_Germine = cor.test(as.numeric(MM2_faces$mm2_z.x),as.numeric(MM2_faces$mm2_z.y))
r_tt_scenes_Germine = cor.test(as.numeric(MM2_scenes$mm2_z.x),as.numeric(MM2_scenes$mm2_z.y))
r_tt_abstracts_Germine = cor.test(as.numeric(MM2_objects$mm2_z.x),as.numeric(MM2_objects$mm2_z.y))

#calculate test-retest icc 2,k (reliability using the average of the two)
icc_tt_faces_Germine = psych::ICC(MM2_faces %>% mutate(mm2_z1 = as.numeric(mm2_z.x), mm2_z2 = as.numeric(mm2_z.y))%>% select(mm2_z1,mm2_z2))
icc_tt_scenes_Germine = psych::ICC(MM2_scenes %>% mutate(mm2_z1 = as.numeric(mm2_z.x), mm2_z2 = as.numeric(mm2_z.y))%>% select(mm2_z1,mm2_z2))
icc_tt_abstracts_Germine = psych::ICC(MM2_objects %>% mutate(mm2_z1 = as.numeric(mm2_z.x), mm2_z2 = as.numeric(mm2_z.y))%>% select(mm2_z1,mm2_z2))

#VALIDATION#####
MM2_TwinLong_val = as.data.frame(MM2_TwinLong_val)

#Prepare MM2_data for MM2 function
names(MM2_TwinLong_val) = c("Obj","Rating","Sub","Block","domain") #rename first 4 columns

MM2_TwinWide_val = MM2_TwinLong_val%>%
  pivot_wider(names_from = "Block", values_from = "Rating")%>%
  filter(!is.na(`2`))

Fac_MM2_TwinWide_val = as.data.frame(MM2_TwinWide_val%>%filter(domain == "FA")%>%
                                   select(-domain))
Sce_MM2_TwinWide_val = as.data.frame(MM2_TwinWide_val%>%filter(domain == "SC")%>%
                                   select(-domain)) 
#extract info from the df
Fac_N_val = length(unique(Fac_MM2_TwinWide_val[,2])) #N participants
Fac_N_items_val = length(unique(Fac_MM2_TwinWide_val[,1]))

#extract info from the df
Sce_N_val = length(unique(Sce_MM2_TwinWide_val[,2])) #N participants
Sce_N_items_val = length(unique(Sce_MM2_TwinWide_val[,1]))

#prepare final matricies for images rated the first time
Fac_MM2_TwinMatrix1_val  = Fac_MM2_TwinWide_val %>%
  select(-`2`) %>% 
  pivot_wider(names_from = "Sub", values_from = "1")
Sce_MM2_TwinMatrix1_val  = Sce_MM2_TwinWide_val %>%
  select(-`2`) %>% 
  pivot_wider(names_from = "Sub", values_from = "1")

####MM2####
#caclulate MM2 (Mean Muinus 2, see MM2.R function)
MM2_faces1_val = MM2(Fac_MM2_TwinMatrix1_val); MM2_faces1_val$MM2
MM2_scenes1_val = MM2(Sce_MM2_TwinMatrix1_val); MM2_scenes1_val$MM2

#prepare final matricies for images rated the second time
Fac_MM2_TwinMatrix2_val  = Fac_MM2_TwinWide_val %>%
  select(-`1`) %>% 
  pivot_wider(names_from = "Sub", values_from = "2")
Sce_MM2_TwinMatrix2_val  = Sce_MM2_TwinWide_val %>%
  select(-`1`) %>% 
  pivot_wider(names_from = "Sub", values_from = "2")

####MM2####
#caclulate MM2 (Mean Muinus 2, see MM2.R function)
MM2_faces2_val = MM2(Fac_MM2_TwinMatrix2_val); MM2_faces2_val$MM2
MM2_scenes2_val = MM2(Sce_MM2_TwinMatrix2_val); MM2_scenes2_val$MM2

#merge mm2 stratified by repeated image
MM2_faces_val = merge(MM2_faces1_val$summary,MM2_faces2_val$summary, by = "Sub")
MM2_scenes_val = merge(MM2_scenes1_val$summary,MM2_scenes2_val$summary, by = "Sub")

#compute test retest on twin one only
MM2_faces_val = MM2_faces_val %>% as.data.frame() %>% 
  mutate(SibID = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibID == "A") %>% 
  select(-SibID)
MM2_scenes_val = MM2_scenes_val %>% as.data.frame() %>% 
  mutate(SibID = substr(Sub,nchar(Sub),nchar(Sub))) %>% 
  filter(SibID == "A") %>% 
  select(-SibID)

#calculate test-retest
r_tt_faces_Sutherland = cor.test(as.numeric(MM2_faces_val$mm2_z.x),as.numeric(MM2_faces_val$mm2_z.y))
r_tt_scenes_Sutherland = cor.test(as.numeric(MM2_scenes_val$mm2_z.x),as.numeric(MM2_scenes_val$mm2_z.y))

#calculate test-retest icc 2,k (reliability using the average of the two)
icc_tt_faces_Sutherland = psych::ICC(MM2_faces_val %>% mutate(mm2_z1 = as.numeric(mm2_z.x), mm2_z2 = as.numeric(mm2_z.y))%>% select(mm2_z1,mm2_z2))
icc_tt_scenes_Sutherland = psych::ICC(MM2_scenes_val %>% mutate(mm2_z1 = as.numeric(mm2_z.x), mm2_z2 = as.numeric(mm2_z.y))%>% select(mm2_z1,mm2_z2))

#SAVE####
save(r_tt_faces_Germine, r_tt_scenes_Germine, r_tt_abstracts_Germine,
     icc_tt_faces_Germine,icc_tt_scenes_Germine, icc_tt_abstracts_Germine,
     r_tt_faces_Sutherland,r_tt_scenes_Sutherland,
     icc_tt_faces_Sutherland,icc_tt_scenes_Sutherland, file = sprintf("%s/%s/04_review/02_report_withinday_reliability_tt.RData",wdOA,wdOA_output))#save output
