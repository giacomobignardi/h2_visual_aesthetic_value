#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#Description: repeate procedure in highlighted in 01_Germine_2015/o4_taste_typicality.R
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
MM2_TwinLong  = read_csv(sprintf("%s/%s/02_Sutherland_2020/01_MM2_Sutherland2020.csv", wdOA,wdOA_output))
#load functions:
source(sprintf("%s/%s/functions/MM2.R", wdOA,wdOA_scripts))

MM2_TwinLong = as.data.frame(MM2_TwinLong)
#Prepare MM2_data for MM2 function
names(MM2_TwinLong) = c("Obj","Rating","Sub","Block","domain") #rename first 4 columns

MM2_TwinWide = MM2_TwinLong%>%
  pivot_wider(names_from = "Block", values_from = "Rating")%>%
  mutate(Rating =  ifelse(is.na(`2`),`1`,(`1`+`2`)/2))%>% #average repeated measures
  select(-c(`1`,`2`))

Fac_MM2_TwinWide = as.data.frame(MM2_TwinWide%>%filter(domain == "FA")%>%
                                   select(-domain))
Sce_MM2_TwinWide = as.data.frame(MM2_TwinWide%>%filter(domain == "SC")%>%
                                   select(-domain)) 
Con_MM2_TwinWide = as.data.frame(MM2_TwinWide%>%filter(domain == "DO")%>%
                                   select(-domain)) 

#extract info from the df
Fac_N = length(unique(Fac_MM2_TwinWide[,2])) #N participants
Fac_N_items = length(unique(Fac_MM2_TwinWide[,1]))

#extract info from the df
Sce_N = length(unique(Sce_MM2_TwinWide[,2])) #N participants
Sce_N_items = length(unique(Sce_MM2_TwinWide[,1]))

#extract info from the df
Con_N = length(unique(Con_MM2_TwinWide[,2])) #N participants
Con_N_items = length(unique(Con_MM2_TwinWide[,1]))

#prepare final matricies
Fac_MM2_TwinMatrix  = Fac_MM2_TwinWide%>%
  pivot_wider(names_from = "Sub", values_from = "Rating")
Sce_MM2_TwinMatrix  = Sce_MM2_TwinWide%>%
  pivot_wider(names_from = "Sub", values_from = "Rating")
Con_MM2_TwinMatrix  = Con_MM2_TwinWide%>%
  pivot_wider(names_from = "Sub", values_from = "Rating")


####MM2####
#caclulate MM2 (Mean Muinus 2, see MM2.R function)
MM2_faces = MM2(Fac_MM2_TwinMatrix); MM2_faces$MM2
MM2_scenes = MM2(Sce_MM2_TwinMatrix); MM2_scenes$MM2
MM2_control = MM2(Con_MM2_TwinMatrix); MM2_control$MM2

####SAVE: dataframe for mm2####
write_csv(MM2_faces$summary%>%mutate(
  mm2_r = as.numeric(mm2_r),
  mm2_z = as.numeric(mm2_z) ),sprintf("%s/%s/02_Sutherland_2020/04_taste_typicality_faces_Sutherland2020.csv",wdOA,wdOA_output))
write_csv(MM2_scenes$summary%>%mutate(
  mm2_r = as.numeric(mm2_r),
  mm2_z = as.numeric(mm2_z) ),sprintf("%s/%s/02_Sutherland_2020/04_taste_typicality_scenes_Sutherland2020.csv",wdOA,wdOA_output))
write_csv(MM2_control$summary%>%mutate(
  mm2_r = as.numeric(mm2_r),
  mm2_z = as.numeric(mm2_z) ),sprintf("%s/%s/02_Sutherland_2020/04_taste_typicality_control_Sutherland2020.csv",wdOA,wdOA_output))