#Author: Giacomo Bignardi
# Last modified: 10-10-2023
#
#
#
#Description:
#collate CTD multivariate results in two comprhensive outputs: 
#1 Saturated Model comparison; 
#2 ACDE models
#clean working environment 
rm(list = ls())
library(tidyverse)
library(psych)
library(readr)
library(patchwork)

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

#load dataFrames:
#multivariate
output_Ger_mul_eb_h  = read_csv(sprintf("%s/%s/01_Germine_2015/13_1_AE_eb_multivariate_h2.csv", wdOA,wdOA_output))
output_Ger_mul_eb_r  = read_csv(sprintf("%s/%s/01_Germine_2015/13_1_AE_eb_multivariate_rArE.csv", wdOA,wdOA_output))

output_Ger_mul_tt_h  = read_csv(sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_h2.csv", wdOA,wdOA_output))
output_Ger_mul_tt_r  = read_csv(sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_rArE.csv", wdOA,wdOA_output))

output_Sut_mul_eb_h  = read_csv(sprintf("%s/%s/02_Sutherland_2020/20_AE_eb_multivariate_Sutherland2020_h2.csv", wdOA,wdOA_output))
output_Sut_mul_eb_r  = read_csv(sprintf("%s/%s/02_Sutherland_2020/20_AE_eb_multivariate_Sutherland2020_rArE.csv", wdOA,wdOA_output))
output_Sut_mul_ebr_h  = read_csv(sprintf("%s/%s/02_Sutherland_2020/22_AE_ebres_multivariate_Sutherland2020_h2.csv", wdOA,wdOA_output))
output_Sut_mul_ebr_r  = read_csv(sprintf("%s/%s/02_Sutherland_2020/22_AE_ebres_multivariate_Sutherland2020_rArE.csv", wdOA,wdOA_output))

output_Sut_mul_tt_h  = read_csv(sprintf("%s/%s/02_Sutherland_2020/21_AE_tt_multivariate_Sutherland2020_h2.csv", wdOA,wdOA_output))
output_Sut_mul_tt_r  = read_csv(sprintf("%s/%s/02_Sutherland_2020/21_AE_tt_multivariate_Sutherland2020_rArE.csv", wdOA,wdOA_output))
output_Sut_mul_ttr_h  = read_csv(sprintf("%s/%s/02_Sutherland_2020/23_AE_ttres_multivariate_Sutherland2020_h2.csv", wdOA,wdOA_output))
output_Sut_mul_ttr_r  = read_csv(sprintf("%s/%s/02_Sutherland_2020/23_AE_ttres_multivariate_Sutherland2020_rArE.csv", wdOA,wdOA_output))

load(sprintf("%s/%s/02_Sutherland_2020/15_sup_FS7_phenotypic_cor_sensitivity_val.rData",wdOA,wdOA_output))
# load(sprintf("%s/%s/02_Sutherland_2020/05_sup_FS7_residuals_sutherland2020.Rdata",wdOA,wdOA_output))
# load(sprintf("%s/%s/02_Sutherland_2020/15_Figure3_phenotypic_cor.rData",wdNOA,wdOA_ImageInput))
# load(sprintf("%s/%s/01_Germine_2015/11_Figure3_phenotypic_cor.rData",wdNOA,wdOA_ImageInput))
# load(sprintf("%s/%s/02_Sutherland_2020/05_supplementary_figureX_residuals_sutherland2020.Rdata",wdNOA,wdOA_ImageInput))
# load(sprintf("%s/%s/02_Sutherland_2020/15_Figure3_phenotypic_cor_sensitivity.rData",wdNOA,wdOA_ImageInput))

# round(output_Ger_mul_eb_h$SA_abstracts[3],2)

####Vizualize####
#Germine####
#Germine heritability estimates
names_germine_x = c("eb abstracts","eb faces","tt scenes","tt abstracts","tt faces")
names_germine_y = c("eb scenes","eb abstracts","eb faces","tt scenes","tt abstracts")
#each SA SE contains h2 or e2 on the diagonal and bh2 or be2 on the off diagonal
output_Ger_eb_SA = output_Ger_mul_eb_h %>% 
  select(SA_scenes,SA_faces,SA_abstracts) %>% 
  mutate(rowname = c("SA_scenes","SA_faces","SA_abstracts"))

h2_eb_abstracts = round(output_Ger_eb_SA %>% filter(rowname == "SA_abstracts") %>% pull(SA_abstracts),2)
h2_eb_scenes = round(output_Ger_eb_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_scenes),2)
h2_eb_faces = round(output_Ger_eb_SA %>% filter(rowname == "SA_faces") %>% pull(SA_faces),2)

bih2_eb_scenesfaces = round(output_Ger_eb_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_faces),2)
bih2_eb_abstractsscenes = round(output_Ger_eb_SA %>% filter(rowname == "SA_abstracts") %>% pull(SA_scenes),2)
bih2_eb_abstractsfaces = round(output_Ger_eb_SA %>% filter(rowname == "SA_abstracts") %>% pull(SA_faces),2)


output_Ger_tt_SA = output_Ger_mul_tt_h %>% 
  select(SA_scenes,SA_faces,SA_abstracts) %>% 
  mutate(rowname = c("SA_scenes","SA_faces","SA_abstracts"))


h2_tt_abstracts = round(output_Ger_tt_SA %>% filter(rowname == "SA_abstracts") %>% pull(SA_abstracts),2)
h2_tt_scenes = round(output_Ger_tt_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_scenes),2)
h2_t_faces = round(output_Ger_tt_SA %>% filter(rowname == "SA_faces") %>% pull(SA_faces),2)

bih2_tt_scenesfaces = round(output_Ger_tt_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_faces),2)

# average_germine_bi = (bih2_eb_scenesfaces + bih2_tt_scenesfaces) / 2

# includes heritability on the diagonal
# Germine = rbind(
#   c(h2_eb_scenes, NA, NA, NA, NA,NA),
#   c(bih2_eb_abstractsscenes, h2_eb_abstracts, NA, NA, NA,NA),
#   c(bih2_eb_scenesfaces, bih2_eb_abstractsfaces, h2_eb_faces, NA, NA,NA),
#   c(0, 0, 0, h2_tt_scenes, NA,NA),
#   c(0, 0, 0, 0, h2_tt_abstracts,NA),
#   c(0, 0, 0, bih2_tt_scenesfaces, 0,h2_tt_faces)
# )

Germine = rbind(
  c(bih2_eb_abstractsscenes, NA, NA, NA,NA),
  c(bih2_eb_scenesfaces, bih2_eb_abstractsfaces, NA, NA,NA),
  c(NA, NA, NA, NA,NA),
  c(NA, NA, NA, NA,NA),
  c(NA, NA, NA, bih2_tt_scenesfaces, NA)
)
colnames(Germine) = names_germine_y
rownames(Germine) = names_germine_x
melted_Germine = reshape2::melt(Germine)

p1_germine = ggplot(data = melted_Germine, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    direction = 1,
    limits = c(0, 1), 
    na.value = 'white')+
  theme_bw(base_size = 18)+
  
  labs(fill = expression(h[b]^2))+  
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  coord_fixed(ratio = 1)

# pdf(sprintf("%s/%s/Figure3b_multivariate_genetic_Germine2015.pdf",wdNOA,wdNOA_ImageOutput),
#     width = 6.5/(1.50),
#     height = 4.5/(1.50))
# p1
# dev.off()

#Germine correlations estimates
#rA
germine_rA_eb_scenesfaces =  round(output_Ger_mul_eb_r %>% filter(correlation == "rA" & trait == "scenes-faces") %>% pull(estimate),2)
germine_rA_eb_abstractsscenes = round(output_Ger_mul_eb_r %>% filter(correlation == "rA" & trait == "scenes-abstracts") %>% pull(estimate),2)
germine_rA_eb_abstractsfaces = round(output_Ger_mul_eb_r %>% filter(correlation == "rA" & trait == "abstracts-faces") %>% pull(estimate),2)
germine_rA_tt_scenesfaces = round(output_Ger_mul_tt_r %>% filter(correlation == "rA") %>% pull(estimate),2)

Germine_rA = rbind(
  c(germine_rA_eb_abstractsscenes, NA, NA, NA,NA),
  c(germine_rA_eb_scenesfaces, germine_rA_eb_abstractsfaces, NA, NA,NA),
  c(NA, NA, NA, NA,NA),
  c(NA, NA, NA, NA ,NA),
  c(NA, NA, NA, germine_rA_tt_scenesfaces, NA)
)
colnames(Germine_rA) = names_germine_y
rownames(Germine_rA) = names_germine_x
melted_Germine_rA = reshape2::melt(Germine_rA)

p2_germine = ggplot(data = melted_Germine_rA, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1, 1),
    na.value = 'white')+
  theme_bw(base_size = 18)+
  
  labs(fill = "ρA")+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("","","","",""))+
  coord_fixed(ratio = 1)

#rE
germine_rE_eb_scenesfaces =  round(output_Ger_mul_eb_r %>% filter(correlation == "rE" & trait == "scenes-faces") %>% pull(estimate),2)
germine_rE_eb_abstractsscenes = round(output_Ger_mul_eb_r %>% filter(correlation == "rE" & trait == "scenes-abstracts") %>% pull(estimate),2)
germine_rE_eb_abstractsfaces = round(output_Ger_mul_eb_r %>% filter(correlation == "rE" & trait == "faces-abstracts") %>% pull(estimate),2)
germine_rE_tt_scenesfaces = round(output_Ger_mul_tt_r %>% filter(correlation == "rE" & trait == "scenes-faces") %>% pull(estimate),2)
germine_rE_tt_abstractsscenes = round(output_Ger_mul_tt_r %>% filter(correlation == "rE" & trait == "scenes-abstracts") %>% pull(estimate),2)
germine_rE_tt_abstractsfaces = round(output_Ger_mul_tt_r %>% filter(correlation == "rE" & trait == "faces-abstracts") %>% pull(estimate),2)


Germine_rE = rbind(
  c(germine_rE_eb_abstractsscenes, NA, NA, NA,NA),
  c(germine_rE_eb_scenesfaces, germine_rE_eb_abstractsfaces, NA, NA,NA),
  c(NA, NA, NA, NA,NA),
  c(NA, NA, NA, germine_rE_tt_abstractsscenes,NA),
  c(NA, NA, NA, germine_rE_tt_scenesfaces, germine_rE_tt_abstractsfaces)
)
colnames(Germine_rE) = names_germine_y
rownames(Germine_rE) = names_germine_x
melted_Germine_rE = reshape2::melt(Germine_rE)

p3_germine = ggplot(data = melted_Germine_rE, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    palette = "BrBG",
    direction = 1,
    limits = c(-1, 1),
    na.value = 'white')+
  theme_bw(base_size = 18)+
  
  labs(fill = "ρE")+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("","","","",""))+
  coord_fixed(ratio = 1)



#Sutherland####
#names = c("eb_scenes","eb_spec_scenes","eb_faces","eb_spec_faces","o_control","tt_scenes","tt_spec_scenes","tt_faces","tt_spec_faces","t_control")
names_sutherland_x = c("eb faces","tt scenes","tt faces")
names_sutherland_y = c("eb scenes","eb faces","tt scenes")

output_Sut_eb_SA = output_Sut_mul_eb_h %>% 
  select(SA_scenes,SA_faces) %>% 
  mutate(rowname = c("SA_scenes","SA_faces"))

output_Sut_ebr_SA = output_Sut_mul_ebr_h %>% 
  select(SA_scenes,SA_faces) %>% 
  mutate(rowname = c("SA_scenes","SA_faces"))


h2_eb_scenes = round(output_Sut_eb_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_scenes),2)
h2_ebr_scenes = round(output_Sut_ebr_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_scenes),2)

h2_eb_faces = round(output_Sut_eb_SA %>% filter(rowname == "SA_faces") %>% pull(SA_faces),2)
h2_ebr_faces = round(output_Sut_ebr_SA %>% filter(rowname == "SA_faces") %>% pull(SA_faces),2)

bih2_eb_scenesfaces = round(output_Sut_eb_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_faces),2)
bih2_ebr_scenesfaces = round(output_Sut_ebr_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_faces),2)

output_Sut_tt_SA = output_Sut_mul_tt_h %>% 
  select(SA_places,SA_faces) %>% 
  mutate(rowname = c("SA_scenes","SA_faces")) %>% 
  rename(SA_scenes = SA_places)

output_Sut_ttr_SA = output_Sut_mul_ttr_h %>% 
  select(SA_scenes,SA_faces) %>% 
  mutate(rowname = c("SA_scenes","SA_faces"))


h2_tt_scenes = round(output_Sut_tt_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_scenes),2)
h2_ttr_scenes = round(output_Sut_ttr_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_scenes),2)

h2_tt_faces = round(output_Sut_tt_SA %>% filter(rowname == "SA_faces") %>% pull(SA_faces),2)
h2_tr_faces = round(output_Sut_ttr_SA %>% filter(rowname == "SA_faces") %>% pull(SA_faces),2)

bih2_tt_scenesfaces = round(output_Sut_tt_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_faces),2)
bih2_ttr_scenesfaces = round(output_Sut_ttr_SA %>% filter(rowname == "SA_scenes") %>% pull(SA_faces),2)

# average_Sutherland_bi = (bih2_eb_scenesfaces + bih2_tt_scenesfaces) / 2
# average_Sutherland_bir = (bih2_ebr_scenesfaces + bih2_ttr_scenesfaces) / 2

sutherland = rbind(
  c(bih2_eb_scenesfaces, NA, NA),
  c(NA, NA,NA),
  c(NA, NA, bih2_tt_scenesfaces)
)

colnames(sutherland) = names_sutherland_y
rownames(sutherland) = names_sutherland_x
melted_sutherland = reshape2::melt(sutherland)

p1_sutherland = ggplot(data = melted_sutherland, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    direction = 1,
    limits = c(0, 1), 
    na.value = 'white')+
  theme_bw(base_size = 18)+
  
  labs(fill = expression(h[b]^2))+  
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  coord_fixed(ratio = 1)


#Sutherland correlations estimates
#rA
#Sutherland genetic correlations estimates
sutherland_rA_eb_scenesfaces =  round(output_Sut_mul_eb_r %>% filter(correlation == "rA" & trait == "scenes-faces") %>% pull(estimate),2)
sutherland_rA_ebr_scenesfaces= round(output_Sut_mul_ebr_r %>% filter(correlation == "rA" & trait == "scenes-faces") %>% pull(estimate),2)
sutherland_rA_tt_scenesfaces =  round(output_Sut_mul_tt_r %>% filter(correlation == "rA" & trait == "places-faces") %>% pull(estimate),2)
sutherland_rA_ttr_scenesfaces = round(output_Sut_mul_ttr_r %>% filter(correlation == "rA" & trait == "scenes-faces") %>% pull(estimate),2)
# 
# average_Sutherland_r = (sutherland_rA_eb_scenesfaces + sutherland_rA_tt_scenesfaces) / 2
# average_Sutherland_rAr = (sutherland_rA_ebr_scenesfaces + sutherland_rA_ttr_scenesfaces) / 2

Sutherland_rA = rbind(
  c(sutherland_rA_eb_scenesfaces, NA, NA),
  c(NA, NA,NA),
  c(NA, NA, sutherland_rA_tt_scenesfaces)
)

colnames(Sutherland_rA) = names_sutherland_y
rownames(Sutherland_rA) = names_sutherland_x
melted_Sutherland_rA = reshape2::melt(Sutherland_rA)


p2_sutherland = ggplot(data = melted_Sutherland_rA, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1, 1),
    na.value = 'white')+
  theme_bw(base_size = 18)+
  labs(fill = "ρA")+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("","",""))+
  coord_fixed(ratio = 1)

#rE
sutherland_rE_eb_scenesfaces =  round(output_Sut_mul_eb_r %>% filter(correlation == "rE" & trait == "scenes-faces") %>% pull(estimate),2)
sutherland_rE_ebr_scenesfaces= round(output_Sut_mul_ebr_r %>% filter(correlation == "rE" & trait == "scenes-faces") %>% pull(estimate),2)
sutherland_rE_tt_scenesfaces =  round(output_Sut_mul_tt_r %>% filter(correlation == "rE" & trait == "places-faces") %>% pull(estimate),2)
sutherland_rE_ttr_scenesfaces = round(output_Sut_mul_ttr_r %>% filter(correlation == "rE" & trait == "scenes-faces") %>% pull(estimate),2)

Sutherland_rE = rbind(
  c(sutherland_rE_eb_scenesfaces, NA, NA),
  c(NA, NA,NA),
  c(NA, NA, sutherland_rE_tt_scenesfaces)
)

colnames(Sutherland_rE) = names_sutherland_y
rownames(Sutherland_rE) = names_sutherland_x
melted_Sutherland_rE = reshape2::melt(Sutherland_rE)

p3_sutherland = ggplot(data = melted_Sutherland_rE, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 3.5)+
  scale_fill_distiller(
    palette = "BrBG",
    direction = 1,
    limits = c(-1, 1),
    na.value = 'white')+
  theme_bw(base_size = 18)+
  labs(fill = "ρE")+
  theme(
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
  scale_y_discrete(labels= c("","",""))+ 
  coord_fixed(ratio = 1)


#Fig5####
pdf(sprintf("%s/%s/04_F5_multivariate_CTD.pdf",wdOA,wdOA_ImageOutput), 
    width = 10,
    height = 6.5)
(p1_germine|p2_germine|p3_germine)/(p1_sutherland|p2_sutherland|p3_sutherland) + plot_layout(guides = 'collect')
dev.off()

#Sensitivity analysis####
names_sutherland_res_x = c("eb res. faces","tt res. scenes","tt res. faces")
names_sutherland_res_y = c("eb res. scenes","eb res. faces","tt res. scenes")

sutherland_r = rbind(
  c(bih2_ebr_scenesfaces, NA, NA),
  c(NA,NA,NA),
  c(NA, NA, bih2_ttr_scenesfaces)
)
colnames(sutherland_r) = names_sutherland_res_y
rownames(sutherland_r) = names_sutherland_res_x
melted_sutherland_r = reshape2::melt(sutherland_r)

Sutherland_rA_r = rbind(
  c(sutherland_rA_ebr_scenesfaces, NA, NA),
  c(NA, NA,NA),
  c(NA, NA, sutherland_rA_ttr_scenesfaces)
)

colnames(Sutherland_rA_r) = names_sutherland_res_y
rownames(Sutherland_rA_r) = names_sutherland_res_x
melted_Sutherland_rA_r = reshape2::melt(Sutherland_rA_r)
Sutherland_rE_r = rbind(
  c(sutherland_rE_ebr_scenesfaces, NA, NA),
  c(NA, NA,NA),
  c(NA, NA, sutherland_rE_ttr_scenesfaces)
)
colnames(Sutherland_rE_r) = names_sutherland_res_y
rownames(Sutherland_rE_r) = names_sutherland_res_x
melted_Sutherland_rE_r = reshape2::melt(Sutherland_rE_r)

#plot bivariate h2; genetic correlations; and environmental correlations for residualised scores
p1_sutherland_r = ggplot(data = melted_sutherland_r, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    direction = 1,
    limits = c(0, 1), 
    na.value = 'white')+
  theme_bw(base_size = 18)+
  labs(fill = expression(h[b]^2))+  
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
   coord_fixed(ratio = 1)

p2_sutherland_r = ggplot(data = melted_Sutherland_rA_r, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 4)+
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1, 1),
    na.value = 'white')+
  theme_bw(base_size = 18)+
  labs(fill = "ρA")+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
   coord_fixed(ratio = 1)
p3_sutherland_r = ggplot(data = melted_Sutherland_rE_r, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white",
            lwd = .2,
            linetype = 1) +
  geom_text(aes(label = round(value,3)), color = "black", size = 3.5)+
  scale_fill_distiller(
    palette = "BrBG",
    direction = 1,
    limits = c(-1, 1),
    na.value = 'white')+
  theme_bw(base_size = 18)+
  labs(fill = "ρE")+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border=element_blank()
  )+
    coord_fixed(ratio = 1)

#FS7####
pdf(sprintf("%s/%s/supplementary/04_SF7_phenotypic_cor_res.pdf",wdOA,wdOA_ImageOutput), 
    width = 9,
    height = 5)
(pCor_t1_sutherland_sensitivity_mm2|pCor_t1_sutherland_sensitivity_eb)+plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a") &theme(legend.position = "right")
dev.off()

#FS8####
pdf(sprintf("%s/%s/supplementary/04_SF8_multivariate_CTD_res.pdf",wdOA,wdOA_ImageOutput), 
    width = 14,
    height = 7)
(p1_sutherland_r|p2_sutherland_r|p3_sutherland_r)+plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a") & theme(legend.position = "none")&theme(legend.position = "right")
dev.off()




