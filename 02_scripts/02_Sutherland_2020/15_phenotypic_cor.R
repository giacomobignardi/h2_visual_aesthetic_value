#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
#compute phenotypic correlations (for twin one and two separately) between evaluation-bias and taste-typicality across domains
#Program: twinDfs ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)
library(umx)
#ggstatsplot

#clean working enviroment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#load dataFrames:
load(file = sprintf("%s/%s/02_Sutherland_2020/05_phenotypic_cor.Rdata",wdOA,wdOA_output))

#Phenotypic correlations
Corr_Pheno = merge(BiometricThin_Pheno_mm2_res_SexAge,BiometricThin_Pheno_avgP_res_SexAge%>%select(-c(AgeAtTestingCalculated, Sex, TwinType)), by = c("Sub","twinN"), all = T)
Corr_Pheno
Corr_Pheno = Corr_Pheno %>%
  rename(
    P_dominance = "DO.y",
    P_faces = "FA.y",
    P_places="SC.y",
    mm2_dominance = "DO.x",
    mm2_faces = "FA.x",
    mm2_places = "SC.x"
  )

Corr_Pheno_mul = merge(BiometricThin_Pheno_mm2_res_multive,BiometricThin_Pheno_avgP_res_multive%>%select(-c(AgeAtTestingCalculated, Sex, TwinType)), by = c("Sub","twinN"), all = T)
Corr_Pheno_mul
Corr_Pheno_mul = Corr_Pheno_mul %>%
  rename(
    Pr_dominance = "DO.y",
    Pr_faces = "FA.y",
    Pr_places="SC.y",
    mm2r_dominance = "DO.x",
    mm2r_faces = "FA.x",
    mm2r_places = "SC.x"
  )

corplot = merge(Corr_Pheno,Corr_Pheno_mul%>%select(-c(AgeAtTestingCalculated, Sex, TwinType)), by = c("Sub","twinN"), all = T)

corplot_t1 = corplot%>%
  filter(twinN == "A")%>%
  select(-c(AgeAtTestingCalculated,Sex, TwinType, twinN, Sub))%>%
  select("P_places","P_faces",
         "mm2_places","mm2_faces")
corplot_t1_sensitivity_p = corplot%>%
  filter(twinN == "A")%>%
  select(-c(AgeAtTestingCalculated,Sex, TwinType, twinN, Sub))%>%
  select("P_places","Pr_places","P_faces","Pr_faces","P_dominance")
corplot_t1_sensitivity_mm2 = corplot%>%
  filter(twinN == "A")%>%
  select(-c(AgeAtTestingCalculated,Sex, TwinType, twinN, Sub))%>%
  select("mm2_places","mm2r_places","mm2_faces","mm2r_faces","mm2_dominance")

#plot correlation for twin one
pCor_t1_sutherland = ggstatsplot::ggcorrmat(
  corplot_t1,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower",
) +
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("eb scenes","eb faces","tt scenes"))+
  scale_x_discrete(labels= c("eb faces","tt scenes","tt faces")) +
  theme(legend.position = "none")


#plot correlation for sensitivity analysis
pCor_t1_sutherland_sensitivity_eb = ggstatsplot::ggcorrmat(
  corplot_t1_sensitivity_p,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower",
) +
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("eb scenes","eb residual\nscenes","eb faces","eb residual\nfaces"))+
  scale_x_discrete(labels= c("eb residual\nscenes","eb faces","eb residual\nfaces", "control"))+
  theme(legend.position = "none")


pCor_t1_sutherland_sensitivity_mm2 = ggstatsplot::ggcorrmat(
  corplot_t1_sensitivity_mm2,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower",
) +
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("tt scenes","tt residual\nscenes","tt faces","tt residual\nfaces"))+
  scale_x_discrete(labels= c("tt residual\nscenes","tt faces","tt residual\nfaces", "control"))+
  theme(legend.position = "none")

corplot_t2 = corplot%>%
  filter(twinN == "B")%>%
  select(-c(AgeAtTestingCalculated,Sex, TwinType, twinN, Sub))%>%
  select("P_places","P_faces",
         "mm2_places","mm2_faces")
corplot_t2_sensitivity_p = corplot%>%
  filter(twinN == "B")%>%
  select(-c(AgeAtTestingCalculated,Sex, TwinType, twinN, Sub))%>%
  select("P_places","Pr_places","P_faces","Pr_faces","P_dominance")
corplot_t2_sensitivity_mm2 = corplot%>%
  filter(twinN == "B")%>%
  select(-c(AgeAtTestingCalculated,Sex, TwinType, twinN, Sub))%>%
  select("mm2_places","mm2r_places","mm2_faces","mm2r_faces","mm2_dominance")

#plot correlation for twin two
pCor_t2_sutherland = ggstatsplot::ggcorrmat(
  corplot_t2,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower",
) +
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("","",""))+
  scale_x_discrete(labels= c("eb faces","tt scenes","tt faces"))+
  theme(legend.position = "none")

#plot correlation for sensitivity analysis
pCor_t2_sutherland_sensitivity_eb = ggstatsplot::ggcorrmat(
  corplot_t2_sensitivity_p,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower",
) +
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("eb scenes","eb residual\nscenes","eb faces","eb residual\nfaces"))+
  scale_x_discrete(labels= c("eb residual\nscenes","eb faces","eb residual\nfaces", "control"))+
  theme(legend.position = "none")


pCor_t2_sutherland_sensitivity_mm2 = ggstatsplot::ggcorrmat(
  corplot_t2_sensitivity_mm2,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower",
) +
  scale_y_discrete(labels= c("tt scenes","tt residual\nscenes","tt faces","tt residual\nfaces"))+
  scale_x_discrete(labels= c("tt residual\nscenes","tt faces","tt residual\nfaces", "control"))+
  theme(legend.position = "none")

#SAVE####
save(pCor_t1_sutherland,pCor_t2_sutherland, file = sprintf("%s/%s/02_Sutherland_2020/15_F3_phenotypic_cor_val.rData",wdOA,wdOA_output))
save(pCor_t1_sutherland_sensitivity_eb,pCor_t1_sutherland_sensitivity_mm2, file = sprintf("%s/%s/02_Sutherland_2020/15_sup_FS7_phenotypic_cor_sensitivity_val.rData",wdOA,wdOA_output))

