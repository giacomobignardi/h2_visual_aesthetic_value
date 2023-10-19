#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description: Compute taste-typicality as the Mean Minus 2 (adapted from Vessel et al.,2018, MM1)
#Program: taste-typicality ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)

#clean working enviroment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#load dataFrames:
MM2_TwinLong  = read_csv(sprintf("%s/%s/01_Germine_2015/01_MM2_Germine2015.csv", wdOA,wdOA_output))
#load functions:
source(sprintf("%s/%s/functions/MM2.R", wdOA,wdOA_scripts))

MM2_TwinLong = as.data.frame(MM2_TwinLong)

#Prepare MM2_data for MM2 function
names(MM2_TwinLong) = c("Obj","Rating","Sub","Block","domain") #rename first 4 columns

MM2_TwinWide = MM2_TwinLong%>%
  pivot_wider(names_from = "Block", values_from = "Rating")%>%
  mutate(Rating =  ifelse(is.na(`2`),`1`,(`1`+`2`)/2))%>% #average repeated measures
  select(-c(`1`,`2`))

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

#prepare final matricies
Fac_MM2_TwinMatrix  = Fac_MM2_TwinWide%>%
  pivot_wider(names_from = "Sub", values_from = "Rating")
Sce_MM2_TwinMatrix  = Sce_MM2_TwinWide%>%
  pivot_wider(names_from = "Sub", values_from = "Rating")
Obj_MM2_TwinMatrix  = Obj_MM2_TwinWide%>%
  pivot_wider(names_from = "Sub", values_from = "Rating")

####MM2####
#caclulate MM2 (Mean Muinus 2, see MM2.R function)
MM2_faces = MM2(Fac_MM2_TwinMatrix); MM2_faces$MM2
MM2_scenes = MM2(Sce_MM2_TwinMatrix); MM2_scenes$MM2
MM2_objects = MM2(Obj_MM2_TwinMatrix); MM2_objects$MM2

####SAVE: dataframe for mm2####
write_csv(MM2_faces$summary%>%mutate(
  mm2_r = as.numeric(mm2_r),
  mm2_z = as.numeric(mm2_z) ),sprintf("%s/%s/01_Germine_2015/04_taste_typicality_faces_Germine2015.csv",wdOA,wdOA_output))
write_csv(MM2_scenes$summary%>%mutate(
  mm2_r = as.numeric(mm2_r),
  mm2_z = as.numeric(mm2_z) ),sprintf("%s/%s/01_Germine_2015/04_taste_typicality_scenes_Germine2015.csv",wdOA,wdOA_output))
write_csv(MM2_objects$summary%>%mutate(
  mm2_r = as.numeric(mm2_r),
  mm2_z = as.numeric(mm2_z) ),sprintf("%s/%s/01_Germine_2015/04_taste_typicality_abstracts_Germine2015.csv",wdOA,wdOA_output))