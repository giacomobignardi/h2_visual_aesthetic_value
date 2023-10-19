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
BioMetric  = read_csv(sprintf("%s/%s/01_Germine_2015/05_twin_multivariate.csv",wdOA,wdOA_output))

#Phenotypic correlations
corplot1_t1 = BioMetric%>%
  select(avgP_1,category,FamId)%>%
  pivot_wider(names_from = "category", values_from = c("avgP_1"))%>%
  rename(P_absObjs = "AO",
         P_faces = "FA_TOT",
         P_scenes = "SC")

corplot2_t1 =BioMetric%>%
  select(mm2_z_1,category,FamId)%>%
  pivot_wider(names_from = "category", values_from = c("mm2_z_1"))%>%
  rename(mm2_absObjs = "AO",
         mm2_faces = "FA_TOT",
         mm2_scenes = "SC")

corplot_t1 = merge(corplot1_t1,corplot2_t1, by = "FamId", all = T)
corplot_t1 = corplot_t1%>%select(P_scenes,P_absObjs,P_faces,mm2_scenes,mm2_absObjs,mm2_faces)

#Plot correlation for twin one
pCor_t1_germine = ggstatsplot::ggcorrmat(
  corplot_t1,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower"
)+
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("eb scenes","eb abstracts","eb faces","tt scenes","tt abstracts"))+
  scale_x_discrete(labels= c("eb abstracts","eb faces","tt scenes","tt abstracts","tt faces"))

#unexpected significant correlation
cor_sce_fac_t1 = cor.test(corplot_t1$mm2_scenes, corplot_t1$P_absObjs)
cor_sce_fac_padj_t1 = p.adjust(cor.test(corplot_t1$mm2_scenes, corplot_t1$P_absObjs)$p.value, "bonferroni",(5*(5+1))/2)
cor_tt_abs_fac_t1 = cor.test(corplot_t1$mm2_faces, corplot_t1$mm2_absObjs)
cor_tt_abs_fac_adj_t1 =p.adjust(cor.test(corplot_t1$mm2_faces, corplot_t1$mm2_absObjs)$p.value, "bonferroni",(5*(5+1))/2)

#####Supplementary Figure 5#####
#repeat on twin second alone
corplot1_t2 = BioMetric%>%
  select(avgP_2,category,FamId)%>%
  pivot_wider(names_from = "category", values_from = c("avgP_2"))%>%
  rename(P_absObjs = "AO",
         P_faces = "FA_TOT",
         P_scenes = "SC")

corplot2_t2 =BioMetric%>%
  select(mm2_z_2,category,FamId)%>%
  pivot_wider(names_from = "category", values_from = c("mm2_z_2"))%>%
  rename(mm2_absObjs = "AO",
         mm2_faces = "FA_TOT",
         mm2_scenes = "SC")

corplot_t2 = merge(corplot1_t2,corplot2_t2, by = "FamId", all = T)
corplot_t2 = corplot_t2%>%select(P_scenes,P_absObjs,P_faces,mm2_scenes,mm2_absObjs,mm2_faces)

#unexpected significant correlation
cor_sce_fac_t2 = cor.test(corplot_t2$mm2_scenes, corplot_t2$P_absObjs)
cor_sce_fac_padj_t2 = p.adjust(cor.test(corplot_t2$mm2_scenes, corplot_t2$P_absObjs)$p.value, "bonferroni",(5*(5+1))/2)
cor_tt_abs_fac_t2 = cor.test(corplot_t2$mm2_faces, corplot_t2$mm2_absObjs)
cor_tt_abs_fac_adj_t2 = p.adjust(cor.test(corplot_t2$mm2_faces, corplot_t2$mm2_absObjs)$p.value, "bonferroni",(5*(5+1))/2)

#save for report
save(cor_sce_fac_t1,cor_sce_fac_padj_t1,
     cor_tt_abs_fac_t1,cor_tt_abs_fac_adj_t1,
     cor_sce_fac_t2,cor_sce_fac_padj_t2,
     cor_tt_abs_fac_t2,cor_tt_abs_fac_adj_t2,
     file = sprintf("%s/%s/01_Germine_2015/11_phenotypic_cor_notsig.rData",wdOA,wdOA_output))


#plot for twin two
pCor_t2_germine = ggstatsplot::ggcorrmat(
  corplot_t2,
  colors = c("#8670b2", "#ffffff", "#f4ab4d"),
  p.adjust.method = "bonferroni",
  type = "Parametric",
  matrix.type = "lower"
)+
  theme_bw(base_size = 18)+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("","","","",""))+
  scale_x_discrete(labels= c("eb abstracts","eb faces","tt scenes","tt abstracts","tt faces"))+
  theme(legend.position = "none")

#SAVE####
save(pCor_t1_germine,pCor_t2_germine, file = sprintf("%s/%s/01_Germine_2015/11_Fig_F3_phenotypic_cor.rData",wdOA,wdOA_output))


