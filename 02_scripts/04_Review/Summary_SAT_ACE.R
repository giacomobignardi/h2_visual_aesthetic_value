#Author: Giacomo Bignardi
#Inspiref from: Verhulst, 2017; Behav Gen
#Date: 2023-06-09
#
#Description: Power analysis for ACE model to detect C component as source of variance in evaluation-bias and taste-typicality
#Program: Power analysis ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(patchwork)
library(tidylog)
library(readr)
require(OpenMx)
require(MASS)

#clean working enviroment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"

#set not Open Access working directories
wdNOA = substr(
  getwd(),
  0,
  nchar(getwd())-nchar("04_analysis_OA")-1
)
wdNOA_Data = "03_rawData/private"
wdNOA_ImageOutput = "05_images/image/processedData"

#SATURATED MODEL AND ASSUMPTION CHECK####
#Discovery####
#Taste-typicality
load(sprintf("%s/%s/01_Germine_2015/08_1_SAT_tt_faces_modelComparison.RData", wdOA,wdOA_output))
modComp_tt_fac = modelComparison
load(sprintf("%s/%s/01_Germine_2015/08_2_SAT_tt_places_modelComparison.RData", wdOA,wdOA_output))
modComp_tt_sce = modelComparison
load(sprintf("%s/%s/01_Germine_2015/08_3_SAT_tt_abstracts_modelComparison.RData", wdOA,wdOA_output))
modComp_tt_abs = modelComparison

#Evaluation-Bias
load(sprintf("%s/%s/01_Germine_2015/07_1_SAT_op_faces_modelComparison.RData", wdOA,wdOA_output))
modComp_eb_fac = modelComparison
load(sprintf("%s/%s/01_Germine_2015/07_2_SAT_op_places_modelComparison.RData", wdOA,wdOA_output))
modComp_eb_sce = modelComparison
load(sprintf("%s/%s/01_Germine_2015/07_3_SAT_op_abstracts_modelComparison.RData", wdOA,wdOA_output))
modComp_eb_abs = modelComparison

#create final df for discovery sample
modComp_dis = rbind(modComp_tt_fac,modComp_tt_sce,modComp_tt_abs,
      modComp_eb_fac,modComp_eb_sce,modComp_eb_abs)

#Evaluation-bias####
#Taste-typicality
load(sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_modelComparison.RData", wdOA,wdOA_output))
modComp_tt_fac_val = modelComparison
load(sprintf("%s/%s/02_Sutherland_2020/08_2_SAT_tt_places_modelComparison.RData", wdOA,wdOA_output))
modComp_tt_sce_val = modelComparison

#Evaluation-Bias
load(sprintf("%s/%s/02_Sutherland_2020/07_1_SAT_op_faces_modelComparison_val.RData", wdOA,wdOA_output))
modComp_eb_fac_val = modelComparison
load(sprintf("%s/%s/02_Sutherland_2020/07_2_SAT_op_places_modelComparison_val.RData", wdOA,wdOA_output))
modComp_eb_sce_val = modelComparison

#create final df for discovery sample
modComp_val = rbind(modComp_tt_fac_val,modComp_tt_sce_val,
                     modComp_eb_fac_val,modComp_eb_sce_val)

#MODEL COMPARISON####
modComp = rbind(modComp_dis,modComp_val)
write_csv(modComp,sprintf("%s/%s/03_review/03_REV_SAT_assumptions.csv",wdOA,wdOA_output))


#FULL ACE MODEL AND MODEL COMPARISON####
#Discovery####
#Taste-typicality
ACE_tt_fac = read_csv(sprintf("%s/%s/01_Germine_2015/10_1_ACE_tt_faces_vc_modelComparsion.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)
ACE_tt_sce = read_csv(sprintf("%s/%s/01_Germine_2015/10_2_ACE_tt_places_vc_modelComparison.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)
ACE_tt_abs = read_csv(sprintf("%s/%s/01_Germine_2015/10_3_ACE_tt_abstracts_vc_modelComparsion.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)


#Evaluation-Bias
ACE_eb_fac = read_csv(sprintf("%s/%s/01_Germine_2015/09_1_ACE_op_faces_vc_modelComparison.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)
ACE_eb_sce = read_csv(sprintf("%s/%s/01_Germine_2015/09_2_ACE_op_places_vc_modelComparison.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)
ADE_eb_abs = read_csv(sprintf("%s/%s/01_Germine_2015/09_3_ADE_op_abstracts_vc_modelComparsion.csv", wdOA,wdOA_output))%>% rename(VCD = VD,SCD = SD)
#create final df for discovery sample
ACE_dis = rbind(ACE_tt_fac,ACE_tt_sce,ACE_tt_abs,
                    ACE_eb_fac,ACE_eb_sce,ADE_eb_abs)

#Evaluation-bias####
#Taste-typicality
ADE_tt_fac_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/12_1_ADE_tt_faces_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VCD = VD,SCD = SD)
ACE_tt_sce_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/12_2_ACE_tt_places_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)

#Evaluation-Bias
ACE_eb_fac_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/11_1_ADE_op_faces_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)
ACE_eb_sce_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/11_2_ADE_op_places_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VCD = VC,SCD = SC)

#create final df for discovery sample
ACE_val = rbind(ADE_tt_fac_val,ACE_tt_sce_val,
                    ACE_eb_fac_val,ACE_eb_sce_val)

#MODEL COMPARISON####
ACE = rbind(ACE_val,ACE_dis)
write_csv(ACE,sprintf("%s/%s/03_review/03_REV_estimates.csv",wdOA,wdOA_output))
