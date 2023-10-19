#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
# Twin Multivariate model to estimate means and (co)variances across multiple groups for residualized data
# Matrix style model - Raw data - Continuous data
#clean working enviroment
rm(list = ls())
library(tidyverse)

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

####DATA####
Fac_eb_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/07_1_SAT_eb_faces_parameters.csv",wdOA,wdOA_output)) %>% mutate(sample = "discovery", domain = "faces", phenotype = "eb")
Sce_eb_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/07_2_SAT_eb_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "discovery", domain = "scenes", phenotype = "eb")
Abs_eb_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/07_3_SAT_eb_abstracts_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "discovery", domain = "abstract images", phenotype = "eb")

Sce_eb_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/07_2_SAT_eb_scenes_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "scenes", phenotype = "eb")
Fac_eb_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/07_1_SAT_eb_faces_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "faces", phenotype = "eb")
Sce_ebr_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/09_2_SAT_ebres_scenes_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "scenes", phenotype = "eb res")
Fac_ebr_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/09_1_SAT_ebres_faces_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "faces", phenotype = "eb res")

Sce_tt_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/08_2_SAT_tt_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "discovery", domain = "scenes", phenotype = "tt")
Fac_tt_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/08_1_SAT_tt_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "discovery", domain = "faces", phenotype = "tt")
Abs_tt_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/08_3_SAT_tt_abstracts_parameters.csv",wdOA,wdOA_output)) %>%mutate(sample = "discovery", domain = "abstract Absects", phenotype = "tt")

Sce_tt_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/08_2_SAT_tt_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "scenes", phenotype = "tt")
Fac_tt_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "faces", phenotype = "tt")
Sce_ttr_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/10_2_SAT_ttres_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "scenes", phenotype = "tt res")
Fac_ttr_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/10_1_SAT_ttres_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "validation", domain = "faces", phenotype = "tt res")

#Table S1####
Twin_phenoR =
rbind(
  #taste typicality
  Abs_tt_Germine,
  Sce_tt_Germine,
  Sce_tt_Sutherland,
  Fac_tt_Germine,
  Fac_tt_Sutherland,
  #evaluation bias
  Abs_eb_Germine,
  Sce_eb_Germine,
  Sce_eb_Sutherland,
  Fac_eb_Germine,
  Fac_eb_Sutherland,
  #taste typicality residuals
  Sce_ttr_Sutherland,
  Fac_ttr_Sutherland,
  #evaluation bias residuals
  Sce_ebr_Sutherland,
  Fac_ebr_Sutherland
      )

view(Twin_phenoR)