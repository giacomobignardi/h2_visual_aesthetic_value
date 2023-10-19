#Author: Giacomo Bignardi
#Date: 18 12 2021
# Last modified: 18-09-2023
#
#
#description: visualize and summarize main output for the article
#NOTE: although this code is labeled 00, it requires some output from higher numbers. Apologies.
#Program: Visualize!!!------------------------------------------------------------------------------------------------------------------------------
#load packages
library(readr)
library(tidyverse)
library(patchwork)

#clear environment
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

####DATA####
#Rxx intra
load(sprintf("%s/%s/01_Germine_2015/01_sup_FS1_exclusion_Germine2015.Rdata",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/01_sup_FS1_exclusion_Sutherland2020.Rdata",wdOA,wdOA_output))

#Phenotypic differences
load(sprintf("%s/%s/01_Germine_2015/06_sup_FS4_tt2pairwise.rData",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/06_sup_FS4_tt2pairwise_val.rData",wdOA,wdOA_output))

#VCA
VCAsummary_main = read_csv(sprintf("%s/%s/01_Germine_2015/02_VCA_twin1_Germine2015.csv",wdOA,wdOA_output))
VCAsummary_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/02_VCA_twin1_Sutherland2020.csv",wdOA,wdOA_output))


VCAsummary_main_t2 = read_csv(sprintf("%s/%s/01_Germine_2015/02_sup_VCA_twin2_Germine2015.csv",wdOA,wdOA_output))
VCAsummary_val_t2 = read_csv(sprintf("%s/%s/02_Sutherland_2020/02_sup_S3_VCA_twin2_Sutherland2020.csv",wdOA,wdOA_output))

#PCA
load(sprintf("%s/%s/01_Germine_2015/05_sup_S5c_PC_Germine2015.Rdata",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/05_sup_S5c_PC_Sutherland2020.Rdata",wdOA,wdOA_output))

load(sprintf("%s/%s/01_Germine_2015/05_sup_FS5_screeplot_PCA.Rdata",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/05_sup_FS5_screeplot_PCA_Sutherland2020.Rdata",wdOA,wdOA_output))

#pairwise agreement
load(sprintf("%s/%s/01_Germine_2015/03_sup_FS6_pairwise_agreement.Rdata",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/03_sup_FS6_pairwise_agreement_val.Rdata",wdOA,wdOA_output))

#Twin correlations
load(sprintf("%s/%s/01_Germine_2015/03_Fig_F2a_pairwise_agreement.Rdata",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/03_Figure2b_pairwise_agreement.Rdata",wdOA,wdOA_output))

Sce_eb_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/07_2_SAT_eb_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "", domain = "scenes", phenotype = "e-b")
Fac_eb_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/07_1_SAT_eb_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "", domain = "faces", phenotype = "e-b") 
Obj_eb_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/07_3_SAT_eb_abstracts_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "", domain = "abstracts", phenotype = "e-b")

Sce_eb_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/07_2_SAT_eb_scenes_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "val.", domain = "scenes", phenotype = "e-b") 
Fac_eb_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/07_1_SAT_eb_faces_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "val.", domain = "faces", phenotype = "e-b") 
Sce_ebres_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/09_2_SAT_ebres_scenes_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "val.r.", domain = "scenes", phenotype = "e-b") 
Fac_ebres_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/09_1_SAT_ebres_faces_parameters_val.csv",wdOA,wdOA_output))%>%mutate(sample = "val.r.", domain = "faces", phenotype = "e-b") 

Sce_mm2_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/08_2_SAT_tt_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "", domain = "scenes", phenotype = "t-t")
Fac_mm2_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/08_1_SAT_tt_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "", domain = "faces", phenotype = "t-t") 
Obj_mm2_Germine = read_csv(sprintf("%s/%s/01_Germine_2015/08_3_SAT_tt_abstracts_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "", domain = "abstracts", phenotype = "t-t")

Sce_mm2_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/08_2_SAT_tt_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "val.", domain = "scenes", phenotype = "t-t")  
Fac_mm2_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "val.", domain = "faces", phenotype = "t-t")  
Sce_mm2res_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/10_2_SAT_ttres_scenes_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "val.r.", domain = "scenes", phenotype = "t-t") 
Fac_mm2res_Sutherland = read_csv(sprintf("%s/%s/02_Sutherland_2020/10_1_SAT_ttres_faces_parameters.csv",wdOA,wdOA_output))%>%mutate(sample = "val.r.", domain = "faces", phenotype = "t-t") 

#h2 comparison
load(sprintf("%s/%s/03_CTD_results/03_CTD_comparison.Rdata",wdOA, wdOA_output))

#phenotypic cross trait correlations
load(sprintf("%s/%s/01_Germine_2015/11_Fig_F3_phenotypic_cor.rData",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/15_F3_phenotypic_cor_val.rData",wdOA,wdOA_output))
load(sprintf("%s/%s/02_Sutherland_2020/15_sup_FS7_phenotypic_cor_sensitivity_val.rData",wdOA,wdOA_output))

# load(sprintf("%s/%s/01_Germine_2015/05_supplementary_figure4_taste_typicality.Rdata",wdOA,wdOA_output))
# load(sprintf("%s/%s/02_Sutherland_2020/05_supplementary_figure4_taste_typicality_val.Rdata",wdOA,wdOA_output))
# load(sprintf("%s/%s/01_Germine_2015/05_supplementary_figure4_overall_pleasantness.Rdata",wdOA,wdOA_output))
# load(sprintf("%s/%s/02_Sutherland_2020/05_supplementary_figure4_overall_pleasantness_val.Rdata",wdOA,wdOA_output))

#multivariate
# output_Ger_mul_eb_h  = read_csv(sprintf("%s/%s/13_1_multivariate_Germine2015_h2.csv", wdOA,wdOA_output))
# output_Ger_mul_eb_r  = read_csv(sprintf("%s/%s/13_1_multivariate_Germine2015_rArE.csv", wdOA,wdOA_output))
# 
# output_Ger_mul_mm2_h  = read_csv(sprintf("%s/%s/13_2_multivariate_Germine2015_h2_mm2.csv", wdOA,wdOA_output))
# output_Ger_mul_mm2_r  = read_csv(sprintf("%s/%s/13_2_multivariate_Germine2015_rArE_mm2.csv", wdOA,wdOA_output))

#F1#####
#rename to plot
VCAsummary_main = VCAsummary_main%>%
  mutate(Component = factor(Component, levels = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock","Residual" )),
         Domain = factor(Domain, levels = c("AO","SC","FA_TOT")), Domain = fct_recode(Domain,"Abstract" = "AO" , "Scenes" = "SC" ,"Faces" = "FA_TOT"))
VCAsummary_val = VCAsummary_val%>%
  mutate(Component = factor(Component, levels = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock","Residual" )),
         Domain = factor(Domain, levels = c("SC","FA")), Domain = fct_recode(Domain, "Scenes val." = "SC" ,"Faces val." = "FA"))

VCAsummary = rbind(VCAsummary_main,VCAsummary_val)

##TEXT: paragraph 1####
#calculate the total amount of individual variance to display in text
VCAsummary %>% filter(grepl("Sub",Component) & !grepl("val",Domain)) %>% group_by(Domain) %>% summarise(sum(Value))

color_VA = c(viridis::viridis(7),c("#D3D3D3"))
#Bar plot with CI
F1b = ggplot(VCAsummary, aes(x=Domain, y=Value, fill = Component)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", alpha = .9) +
  geom_linerange(aes(ymin=CI_low, ymax=CI_high), position=position_dodge(.9),size = 1) +
  scale_fill_manual(values = color_VA) +
  scale_y_continuous(limits = c(0,.70), breaks = seq(0,.7,by = .1))+
  labs(y= "Proportion of \n variance", x = "Domain" ,fill = "Variance\n Component")+
  theme_classic(base_size = 18)

##F1b#####
pdf(sprintf("%s/%s/00_F1b_VPC.pdf",wdOA,wdOA_ImageOutput),
    width = 10,
    height = 3.5)
F1b
dev.off()


#F2####
Twin_phenoR =
  rbind(
    Sce_eb_Germine,Fac_eb_Germine,Obj_eb_Germine,
    Sce_eb_Sutherland,Fac_eb_Sutherland,Sce_ebres_Sutherland,Fac_ebres_Sutherland,
    
    Sce_mm2_Germine,Fac_mm2_Germine,Obj_mm2_Germine,
    Sce_mm2_Sutherland,Fac_mm2_Sutherland,Sce_mm2res_Sutherland,Fac_mm2res_Sutherland
  )

#create a variable to convienently plot correlation in 2d
Twin_pheno = Twin_phenoR%>%
  filter( substr(parameter, 1, 2) == "MZ" | substr(parameter, 1, 2) == "DZ")%>%
  mutate(
    Zygosity = substr(parameter, 1, 2),
    Domain = sprintf("%s %s",domain,sample))

#col = viridis::viridis(3)[2:3] col before review  
col_updt = wesanderson::wes_palette("Darjeeling1")[c(2,4,5)]

F2c  = ggplot(Twin_pheno%>%filter(sample == "")%>%filter(phenotype == "t-t"), aes(domain, y = estimate, color = Zygosity)) + 
  geom_point(size = 2, position = position_dodge(.9)) +
  geom_errorbar(aes(ymin = lbound, ymax = ubound), size = 1, width = .1, position = position_dodge(.9)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.1) +
  scale_color_manual(values = rev(col_updt[1:2])) + 
  ylim(-.1,1) +
  theme_classic(base_size = 14) +
  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5),
        legend.position = "none") +
  scale_x_discrete(labels= c("abstracts","faces","scenes"))+
  labs(y = "",
       x = "Taste-typicality")

F2e  = ggplot(Twin_pheno%>%filter(sample == "")%>%filter(phenotype == "e-b"), aes(domain, y = estimate, color = Zygosity)) + 
  geom_point(size = 2, position = position_dodge(.9)) +
  geom_errorbar(aes(ymin = lbound, ymax = ubound), size = 1, width = .1, position = position_dodge(.9)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.1) +
  scale_color_manual(values = rev(col_updt[1:2])) + 
  ylim(-.1,1) +
  theme_classic(base_size = 14) +
  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5),
        legend.position = "none") +
  scale_x_discrete(labels= c("abstracts","faces","scenes"))+
  labs(y = "Phenotypic r",
       x = "Evaluation-bias")

F2d  = ggplot(Twin_pheno%>%
                filter(sample == "val.")%>%filter(phenotype == "t-t"), aes(Domain, y = estimate, color = Zygosity)) + 
  geom_point(size = 2, position = position_dodge(.9)) +
  geom_errorbar(aes(ymin = lbound, ymax = ubound), size = 1, width = .1, position = position_dodge(.9)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.1) +
  scale_color_manual(values = rev(col_updt[1:2])) + 
  ylim(-.1,1) +
  theme_classic(base_size = 14) +
  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5),
        legend.position = "none") +
  scale_x_discrete(labels= c("faces","scenes"))+
  labs(y = "",
       x = "")

F2f  = ggplot(Twin_pheno%>%filter(sample == "val.")%>%filter(phenotype == "e-b"), aes(Domain, y = estimate, color = Zygosity)) + 
  geom_point(size = 2, position = position_dodge(.9)) +
  geom_errorbar(aes(ymin = lbound, ymax = ubound), size = 1, width = .1, position = position_dodge(.9)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.1) +
  scale_color_manual(values = rev(col_updt[1:2])) + 
  ylim(-.1,1) +
  theme_classic(base_size = 14) +
  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5),
        legend.position = "none") +
  scale_x_discrete(labels= c("faces","scenes"))+
  labs(y = "Phenotypic r",
       x = "")

F2a = p_pairwise_r + scale_color_manual(values = rev(col_updt)) +  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5)) + ylab("pairwise-agreement")+ xlab("")+theme(legend.position = "bottom") + scale_x_discrete(labels= c("abstracts","scenes","faces"))
F2b = p_pairwise_r_val + scale_color_manual(values = rev(col_updt)) +  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5))+ylab("")+xlab("")+  scale_x_discrete(labels= c("scenes","faces"))+ theme(legend.position = "bottom") + scale_x_discrete(labels= c("scenes","faces"))


pdf(sprintf("%s/%s/00_F2_CTD.pdf",wdOA,wdOA_ImageOutput), 
    width = 13.5, 
    height = 8)
((((F2a | F2b) /
    (F2c | F2d) /
    (F2e | F2f))+ plot_layout(guides = 'collect')) |(h2_comparison +theme_classic(base_size = 14) + theme(legend.position = "bottom"))) + plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(face = 'bold'),legend.position = 'bottom')
dev.off()

#F3####
pdf(sprintf("%s/%s/00_F3_crosstrait_pheno_cor.pdf",wdOA,wdOA_ImageOutput), 
    width = 13.5, 
    height = 4)
(pCor_t1_germine|pCor_t2_germine|pCor_t1_sutherland|pCor_t2_sutherland)+ plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
dev.off()

#F4####
#see 04 CTD multivariate for figure 4

#SUPPLEMENTARY####
#SF1:Rxx####
p1_exclusion_abs = p1_exclusion_abs + xlab("Rxx-intra r abstracts (discovery)")
p1_exclusion_pla = p1_exclusion_pla + xlab("Rxx-intra r scenes (discovery)")
p1_exclusion_fac =p1_exclusion_fac + xlab("Rxx-intra r faces (discovery)")
p1_exclusion_pla_val =p1_exclusion_pla_val + xlab("Rxx-intra r scenes (validation)")
p1_exclusion_fac_val =p1_exclusion_fac_val + xlab("Rxx-intra r faces (validation)")

pdf(sprintf("%s/%s/supplementary/SF1_exclusion.pdf",wdOA,wdOA_ImageOutput), 
    width = 9.5, 
    height = 5)
(p1_exclusion_abs + p1_exclusion_pla + p1_exclusion_fac)/(plot_spacer() + p1_exclusion_pla_val + p1_exclusion_fac_val) + plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(face = 'bold'))
dev.off()

#SF3####
##SF3a####
#rename to plot
VCAsummary_main_t2_sup = VCAsummary_main_t2%>%
  mutate(Component = factor(Component, levels = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock","Residual" )),
         Component = fct_recode(Component,Image = "Obj", Individual = "Sub", Exposure = "Block", ImageXindividual = "SubXObj", ExposureXindividual = "SubXBlock",  ImageXExposure = "ObjxBlock", Residual = "Residual"),
         Domain = factor(Domain, levels = c("AO","SC","FA_TOT")),
         Domain = fct_recode(Domain,"Abstract" = "AO" , "Scenes" = "SC" ,"Faces" = "FA_TOT"))
VCAsummary_main_t2_sup %>% group_by(Component,Domain) %>% summarise(Value = mean(Value)) %>% filter(Component == "Individual"| Component == "ImageXindividual") %>% group_by(Domain) %>% summarise(mean(Value))

VCAsummary_val_t2_sup = VCAsummary_val_t2%>%
  mutate(Component = factor(Component, levels = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock","Residual" )),
         Component = fct_recode(Component,Image = "Obj", Individual = "Sub", Exposure = "Block", ImageXindividual = "SubXObj", ExposureXindividual = "SubXBlock",  ImageXExposure = "ObjxBlock", Residual = "Residual"),
         Domain = factor(Domain, levels = c("SC","FA")), Domain = fct_recode(Domain, "Scenes val." = "SC" ,"Faces val." = "FA"))

VCAsummary_t2_sup = rbind(VCAsummary_main_t2_sup,VCAsummary_val_t2_sup)

color_VA = c(viridis::viridis(7),c("#D3D3D3"))
#Bar plot with CI
SF2a = ggplot(VCAsummary_t2_sup, aes(x=Domain, y=Value, fill = Component)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", alpha = .9) +
  geom_linerange(aes(ymin=CI_low, ymax=CI_high), position=position_dodge(.9),size = 1) +
  scale_fill_manual(values = color_VA) +
  scale_y_continuous(limits = c(0,.70), breaks = seq(0,.7,by = .1))+
  labs(y= "Proportion of \n variance", x = "Domain" ,fill = "Variance\n Component")+
  theme_classic(base_size = 16)
####SF3b####
#vizualize Honnekop modified beolder index
VCAsummary_Bi = VCAsummary%>%
  mutate(Shared = 1-mBi)%>%
  rename(Unique = "mBi")%>%
  dplyr::select(Domain, Unique, Shared)%>%
  pivot_longer(c("Unique", "Shared"), names_to = "Beholder_Index", values_to = "Value")%>%
  group_by(Beholder_Index,Domain)%>%
  summarise(Value = mean(Value))

SF2b = ggplot(VCAsummary_Bi, aes(x= Value, y=Domain, fill=Beholder_Index))+
  geom_col(width = .7) +
  scale_fill_manual(values = c("#440154FF","#316988"), labels = c("1 - mB1 (shared)", "mB1 (unique)")) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "gray", size = 1) +
  labs(
    fill = "Beholder Index",
    x = "Proportion of variance",
    y = " Domain",
    subtitle = "twin 1"
  )+
  theme_classic(base_size = 16)

VCAsummary_Bi_t2 = VCAsummary_t2_sup%>%
  mutate(Shared = 1-mBi)%>%
  rename(Unique = "mBi")%>%
  dplyr::select(Domain, Unique, Shared)%>%
  pivot_longer(c("Unique", "Shared"), names_to = "Beholder_Index", values_to = "Value")%>%
  group_by(Beholder_Index,Domain)%>%
  summarise(Value = mean(Value))

SF2b_t2 = ggplot(VCAsummary_Bi_t2, aes(x= Value, y=Domain, fill=Beholder_Index))+
  geom_col(width = .7) +
  scale_fill_manual(values = c("#440154FF","#316988"), labels = c("1 - mB1 (shared)", "mB1 (unique)")) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "gray", size = 1) +
  labs(
    fill = "Beholder Index",
    x = "Proportion of variance",
    y = " Domain",
    subtitle = "twin 2"
  )+
  theme_classic(base_size = 16)

# #save
pdf(sprintf("%s/%s/supplementary/00_SF3_VCA_t2.pdf",wdOA,wdOA_ImageOutput), 
    width = 11,
    height = 6)
((SF2a) / ((SF2b |SF2b_t2)+
             plot_layout(guides = 'collect'))) + plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(face = 'bold'))
dev.off()

#SF4:phenotypic comparison (tt 2 pa)#####
p_dif_abstract=p_dif_abstract + ylab("pairwise agreement") + xlab(expression(abs(Delta[paste(taste-typicality,sep ="_",abstract)])))
p_dif_scenes=p_dif_scenes + ylab("pairwise agreement") + xlab(expression(abs(Delta[paste(taste-typicality,sep ="_",scenes)])))
p_dif_faces=p_dif_faces + ylab("pairwise agreement") + xlab(expression(abs(Delta[paste(taste-typicality,sep ="_",faces)])))
p_dif_scenes_val=p_dif_scenes_val + ylab("pairwise agreement val.") + xlab(expression(abs(Delta[paste(taste-typicality,sep ="_",scenes,sep =" ",val)])))
p_dif_faces_val=p_dif_faces_val + ylab("pairwise agreement val.") + xlab(expression(abs(Delta[paste(taste-typicality,sep ="_",faces,sep =" ",val)])))

pdf(sprintf("%s/%s/supplementary/SF4_tasteTypicality_pairwiseAgreement.pdf",wdOA,wdOA_ImageOutput), 
    width = 9, 
    height = 6)
(p_dif_abstract + p_dif_scenes + p_dif_faces)/(plot_spacer() + p_dif_scenes_val + p_dif_faces_val) + plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(face = 'bold'))
dev.off()

#SF5:PCA#####
Figure_SF5_1 = ((p_sp_1| p_sp_2 |p_sp_1_sut2020| p_sp_2_sut2020)+ plot_layout(guides = "collect"))

Figure_SF5_2 = ((p_PC2_AO|p_PC1_AO)/(p_PC2_SC|p_PC1_SC)/(p_PC2_FA|p_PC1_FA) | 
  ((plot_spacer()|plot_spacer())/(p_PC2_SC_sut2020|p_PC1_SC_sut2020)/(p_PC2_FA_sut2020|p_PC1_FA_sut2020)) ) + plot_layout(guides = "keep") 
Figure_SF5_3 = (p_Obj_BioMetric_PC_avg_t2_Ger2015 | p_Sce_BioMetric_PC_avg_t2_Ger2015 | p_Fac_BioMetric_PC_avg_t2_Ger2015 |p_Sce_BioMetric_PC_avg_t2_Sut2020 | p_Fac_BioMetric_PC_avg_t2_Sut2020)/
  (p_Obj_BioMetric_PC_mm2_t2_Ger2015 | p_Sce_BioMetric_PC_mm2_t2_Ger2015 | p_Fac_BioMetric_PC_mm2_t2_Ger2015 |p_Sce_BioMetric_PC_mm2_t2_Sut2020 | p_Fac_BioMetric_PC_mm2_t2_Sut2020)

pdf(sprintf("%s/%s/supplementary/SF5a_PCA.pdf",wdOA,wdOA_ImageOutput), 
    width = 9,
    height = 2.5)
Figure_SF5_1
dev.off()
pdf(sprintf("%s/%s/supplementary/SF5b_PCA.pdf",wdOA,wdOA_ImageOutput), 
    width = 20,
    height = 10)
Figure_SF5_2
dev.off()
pdf(sprintf("%s/%s/supplementary/SF5c_PCA.pdf",wdOA,wdOA_ImageOutput), 
    width = 11.5,
    height = 4)
Figure_SF5_3
dev.off()

#SF6: pairwise agreement####
p_pairwise_z= p_pairwise_z + 
  scale_color_manual(values = rev(col_updt))+
  theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5)) +
  ylab("") +
  scale_x_discrete(labels= c("abstracts","scenes","faces"))
p_pairwise_z_val= p_pairwise_z_val + theme(axis.text.x =  element_text(angle = 45, hjust = .5, vjust = .5))+ylab("")+  scale_color_manual(values = rev(col_updt))
p1= p1  + ylab("pa") + scale_fill_manual(values = rev(col_updt))
p2= p2  + ylab("pa (Fisher-z)") +
  scale_x_continuous(
    breaks = c(0,1,2,3)
  ) + scale_fill_manual(values = rev(col_updt))

p1_val= p1_val + 
  ylab("pa") +
  scale_x_continuous(
  breaks = scales::pretty_breaks(),
)+ scale_fill_manual(values = rev(col_updt))
p2_val= p2_val  + 
  ylab("pa (Fisher-z)") +
  scale_x_continuous(
    breaks = c(0,1,2,3)
  )+ scale_fill_manual(values = rev(col_updt))

FS6_1 = ((((p1|p2) +plot_layout(guides = 'collect'))|p_pairwise_z)/
    (p3|p4|p5)) + plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(face = 'bold'))

FS6_2 = ((((p1_val|p2_val) +plot_layout(guides = 'collect'))|p_pairwise_z_val)/
    (p3_val|p4_val|p5_val)) + plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(face = 'bold'))

pdf(sprintf("%s/%s/supplementary/SF6a_pairwise_agreement.pdf",wdOA,wdOA_ImageOutput), 
    width = 12, 
    height = 11)
FS6_1
dev.off()

pdf(sprintf("%s/%s/supplementary/SF6b_pairwise_agreement.pdf",wdOA,wdOA_ImageOutput), 
    width = 12, 
    height = 11)
FS6_2
dev.off()