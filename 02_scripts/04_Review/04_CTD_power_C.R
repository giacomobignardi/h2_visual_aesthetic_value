#Author: Giacomo Bignardi
#Inspired by: Verhulst, 2017; Behav Gen
#Date: 2023-06-09
#Last modified: 2023-10-10
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

#clean working environment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

#load functions:
source(sprintf("%s/%s/functions/powerFun.R", wdOA,wdOA_scripts)) #Verhulst, 2017; Behav Gen https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471839/
source(sprintf("%s/%s/functions/ACE_powCurve.R", wdOA,wdOA_scripts))

output = read_csv(sprintf("%s/%s/03_CTD_results/02_SF2_ACE_results.csv", wdOA,wdOA_output))
#select A and C components
output_C = output %>% dplyr::select(domain, facet, comparison, sample,SC_SD) %>% filter(comparison == "ACE")
output_A = output %>% dplyr::select(domain, facet, comparison, sample,SA) %>% filter(comparison == "ACE")

#Power-analysis: GENERATED DATA####
#Generated data: Step follow procedure as in Verhulst, 2017; Behav Gen
#Overiew of steps:
##Step 1: simulate twin data
##Step 2: fit full and reduced model and obtain the χ2 for the LRT
##Step 3: calculate the Weighted Non Centrality Parameter χ2/(nMZp + nDZp), where p stands for pair
##step 1 to 3 is achieved in the Pow function in powerFun: here is based on mean twin-h2 estimated for each trait, across different values of C

#Number of twin based on Table 1
nMZ_dis = 133 + 424  
nDZ_dis = 54 + 162

nMZ_val = 94 + 245  
nDZ_val = 23 + 137

#total family-wise sample
x_sample_dis = nMZ_dis + nDZ_dis
x_sample_val = nMZ_val + nDZ_val

#ACE model
output_C %>% dplyr::select(domain, facet, sample)
#Step 1 to 3
#extract A estimates
quant_dis_sce_eb = c()
A_dis_sce_eb  = output_A %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Germine et al. 2015" ) %>% pull(SA)
quant_dis_fac_eb = c()
A_dis_fac_eb  = output_A %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Germine et al. 2015" ) %>% pull(SA)

quant_dis_abs_tt = c()
A_dis_abs_tt  = output_A %>% filter(domain == "abstract" & facet == "taste_typicality", sample == "Germine et al. 2015" ) %>% pull(SA)
quant_dis_sce_tt = c()
A_dis_sce_tt  = output_A %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Germine et al. 2015" ) %>% pull(SA)
quant_dis_fac_tt = c()
A_dis_fac_tt  = output_A %>% filter(domain == "faces" & facet == "taste_typicality", sample == "Germine et al. 2015" ) %>% pull(SA)

quant_val_sce_eb = c()
A_val_sce_eb  = output_A %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" ) %>% pull(SA)
quant_val_fac_eb = c()
A_val_fac_eb  = output_A %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" ) %>% pull(SA)

quant_val_sce_tt = c()
A_val_sce_tt  = output_A %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Sutherland et al. 2020" ) %>% pull(SA)

#create a sequence of C to generate several power curves
C = seq(.01,.35,.01)
#generate simulated WNCP for each combination of A(esitmated,fixed) and C(.01 to .35)
for(i in C){
quant_dis_sce_eb_i = acePow( A_dis_sce_eb, i, nMZ_dis, nDZ_dis)
quant_dis_fac_eb_i = acePow( A_dis_fac_eb, i, nMZ_dis, nDZ_dis)

quant_dis_abs_tt_i = acePow( A_dis_abs_tt, i, nMZ_dis, nDZ_dis)
quant_dis_sce_tt_i = acePow( A_dis_sce_tt, i, nMZ_dis, nDZ_dis)
quant_dis_fac_tt_i = acePow( A_dis_fac_tt, i, nMZ_dis, nDZ_dis)

quant_val_sce_eb_i = acePow( A_val_sce_eb, i, nMZ_val, nDZ_val)
quant_val_fac_eb_i = acePow( A_val_fac_eb, i, nMZ_val, nDZ_val)

quant_val_sce_tt_i = acePow( A_val_sce_tt, i, nMZ_val, nDZ_val)

quant_dis_sce_eb = rbind(quant_dis_sce_eb, data.frame(WncpA = quant_dis_sce_eb_i$WncpA,WncpC = quant_dis_sce_eb_i$WncpC, C = i))
quant_dis_fac_eb = rbind(quant_dis_fac_eb, data.frame(WncpA = quant_dis_fac_eb_i$WncpA,WncpC = quant_dis_fac_eb_i$WncpC, C = i))

quant_dis_abs_tt = rbind(quant_dis_abs_tt, data.frame(WncpA = quant_dis_abs_tt_i$WncpA,WncpC = quant_dis_abs_tt_i$WncpC, C = i))
quant_dis_sce_tt = rbind(quant_dis_sce_tt, data.frame(WncpA = quant_dis_sce_tt_i$WncpA,WncpC = quant_dis_sce_tt_i$WncpC, C = i))
quant_dis_fac_tt = rbind(quant_dis_fac_tt, data.frame(WncpA = quant_dis_fac_tt_i$WncpA,WncpC = quant_dis_fac_tt_i$WncpC, C = i))

quant_val_sce_eb = rbind(quant_val_sce_eb, data.frame(WncpA = quant_val_sce_eb_i$WncpA,WncpC = quant_val_sce_eb_i$WncpC, C = i))
quant_val_fac_eb = rbind(quant_val_fac_eb, data.frame(WncpA = quant_val_fac_eb_i$WncpA,WncpC = quant_val_fac_eb_i$WncpC, C = i))

quant_val_sce_tt = rbind(quant_val_sce_tt, data.frame(WncpA = quant_val_sce_tt_i$WncpA,WncpC = quant_val_sce_tt_i$WncpC, C = i))
}

#mutate C to factor for later plotting
quant_dis_sce_eb = quant_dis_sce_eb %>% mutate(C = factor(C))
quant_dis_fac_eb = quant_dis_fac_eb %>% mutate(C = factor(C))
quant_dis_abs_tt = quant_dis_abs_tt %>% mutate(C = factor(C))
quant_dis_sce_tt = quant_dis_sce_tt %>% mutate(C = factor(C))
quant_dis_fac_tt = quant_dis_fac_tt %>% mutate(C = factor(C))
quant_val_sce_eb = quant_val_sce_eb %>% mutate(C = factor(C))
quant_val_fac_eb = quant_val_fac_eb %>% mutate(C = factor(C))
quant_val_sce_tt = quant_val_sce_tt %>% mutate(C = factor(C))

powCurve_dis_sce_eb = c()
powCurve_dis_fac_eb = c()
powCurve_dis_abs_tt = c()
powCurve_dis_sce_tt = c()
powCurve_dis_fac_tt = c()
powCurve_val_sce_eb = c()
powCurve_val_fac_eb = c()
powCurve_val_sce_tt = c()

#Step 4: get WNCP per unit of increment and compute power based on the χ2(1) test (based on LRT)
#Multiply WNCP per number of family and compute power
for(i in C){
powCurve_dis_sce_eb_i = cbind(powCurveACE(quant_dis_sce_eb %>% filter(C == i)),C = i)
powCurve_dis_fac_eb_i = cbind(powCurveACE(quant_dis_fac_eb %>% filter(C == i)),C = i)
powCurve_dis_abs_tt_i = cbind(powCurveACE(quant_dis_abs_tt %>% filter(C == i)),C = i)
powCurve_dis_sce_tt_i = cbind(powCurveACE(quant_dis_sce_tt %>% filter(C == i)),C = i)
powCurve_dis_fac_tt_i = cbind(powCurveACE(quant_dis_fac_tt %>% filter(C == i)),C = i)
powCurve_val_sce_eb_i = cbind(powCurveACE(quant_val_sce_eb %>% filter(C == i)),C = i)
powCurve_val_fac_eb_i = cbind(powCurveACE(quant_dis_fac_eb %>% filter(C == i)),C = i)
powCurve_val_sce_tt_i = cbind(powCurveACE(quant_dis_sce_tt %>% filter(C == i)),C = i)
powCurve_dis_sce_eb = rbind(powCurve_dis_sce_eb,powCurve_dis_sce_eb_i)
powCurve_dis_fac_eb = rbind(powCurve_dis_fac_eb,powCurve_dis_fac_eb_i)
powCurve_dis_abs_tt = rbind(powCurve_dis_abs_tt,powCurve_dis_abs_tt_i)
powCurve_dis_sce_tt = rbind(powCurve_dis_sce_tt,powCurve_dis_sce_tt_i)
powCurve_dis_fac_tt = rbind(powCurve_dis_fac_tt,powCurve_dis_fac_tt_i)
powCurve_val_sce_eb = rbind(powCurve_val_sce_eb,powCurve_val_sce_eb_i)
powCurve_val_fac_eb = rbind(powCurve_val_fac_eb,powCurve_val_fac_eb_i)
powCurve_val_sce_tt = rbind(powCurve_val_sce_tt,powCurve_val_sce_tt_i)
}

#Power-analysis:  EMPRICAL WNCP ####
#Step 1: obtain the difference between the ACE and the AE model
X2C_dis_sce_eb =  output %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Germine et al. 2015" , comparison  == "AE") %>% pull(diffLL)
X2A_dis_sce_eb =  output %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Germine et al. 2015" , comparison  == "CE") %>% pull(diffLL)
X2C_dis_fac_eb =  output %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Germine et al. 2015" , comparison  == "AE") %>% pull(diffLL)
X2A_dis_fac_eb =  output %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Germine et al. 2015" , comparison  == "CE") %>% pull(diffLL)

X2C_dis_abs_tt =  output %>% filter(domain == "abstract" & facet == "taste_typicality", sample == "Germine et al. 2015" , comparison  == "AE") %>% pull(diffLL)
X2A_dis_abs_tt =  output %>% filter(domain == "abstract" & facet == "taste_typicality", sample == "Germine et al. 2015" , comparison  == "CE") %>% pull(diffLL)
X2C_dis_sce_tt =  output %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Germine et al. 2015" , comparison  == "AE") %>% pull(diffLL)
X2A_dis_sce_tt =  output %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Germine et al. 2015" , comparison  == "CE") %>% pull(diffLL)
X2C_dis_fac_tt =  output %>% filter(domain == "faces" & facet == "taste_typicality", sample == "Germine et al. 2015" , comparison  == "AE") %>% pull(diffLL)
X2A_dis_fac_tt =  output %>% filter(domain == "faces" & facet == "taste_typicality", sample == "Germine et al. 2015" , comparison  == "CE") %>% pull(diffLL)

X2C_val_sce_eb =  output %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" , comparison  == "AE") %>% pull(diffLL)
X2A_val_sce_eb =  output %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" , comparison  == "CE") %>% pull(diffLL)
X2C_val_fac_eb =  output %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" , comparison  == "AE") %>% pull(diffLL)
X2A_val_fac_eb =  output %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" , comparison  == "CE") %>% pull(diffLL)

X2C_val_sce_tt =  output %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Sutherland et al. 2020" , comparison  == "AE") %>% pull(diffLL)
X2A_val_sce_tt =  output %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Sutherland et al. 2020" , comparison  == "CE") %>% pull(diffLL)



#Step 2: divide by the number of total observation (family wise, for simplicity we use the full sample)
Wncp_dis_sce_eb  = data.frame(WncpC = X2C_dis_sce_eb/x_sample_dis,WncpA = X2A_dis_sce_eb/x_sample_dis)
Wncp_dis_fac_eb  = data.frame(WncpC = X2C_dis_fac_eb/x_sample_dis,WncpA = X2A_dis_fac_eb/x_sample_dis)

Wncp_dis_abs_tt  = data.frame(WncpC = X2C_dis_abs_tt/x_sample_dis,WncpA = X2A_dis_abs_tt/x_sample_dis)
Wncp_dis_sce_tt  = data.frame(WncpC = X2C_dis_sce_tt/x_sample_dis,WncpA = X2A_dis_sce_tt/x_sample_dis)
Wncp_dis_fac_tt  = data.frame(WncpC = X2C_dis_fac_tt/x_sample_dis,WncpA = X2A_dis_fac_tt/x_sample_dis)

Wncp_val_sce_eb  = data.frame(WncpC = X2C_val_sce_eb/x_sample_val,WncpA = X2A_val_sce_eb/x_sample_val)
Wncp_val_fac_eb  = data.frame(WncpC = X2C_val_fac_eb/x_sample_val,WncpA = X2A_val_fac_eb/x_sample_val)

Wncp_val_sce_tt  = data.frame(WncpC = X2C_val_sce_tt/x_sample_val,WncpA = X2A_val_sce_tt/x_sample_val)

#Step 3: obtain value for a range of sample sizes
powCurve_dis_sce_eb_obs = powCurveACE(Wncp_dis_sce_eb)
powCurve_dis_fac_eb_obs = powCurveACE(Wncp_dis_fac_eb)

powCurve_dis_abs_tt_obs = powCurveACE(Wncp_dis_abs_tt)
powCurve_dis_sce_tt_obs = powCurveACE(Wncp_dis_sce_tt)
powCurve_dis_fac_tt_obs = powCurveACE(Wncp_dis_fac_tt)

powCurve_val_sce_eb_obs = powCurveACE(Wncp_val_sce_eb)
powCurve_val_fac_eb_obs = powCurveACE(Wncp_val_fac_eb)

powCurve_val_sce_tt_obs = powCurveACE(Wncp_val_sce_tt)


#extract estimated C
C_dis_sce_eb  = output_C %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Germine et al. 2015" ) %>% pull(SC_SD)
C_dis_fac_eb  = output_C %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Germine et al. 2015" ) %>% pull(SC_SD)

C_dis_abs_tt  = output_C %>% filter(domain == "abstract" & facet == "taste_typicality", sample == "Germine et al. 2015" ) %>% pull(SC_SD)
C_dis_sce_tt  = output_C %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Germine et al. 2015" ) %>% pull(SC_SD)
C_dis_fac_tt  = output_C %>% filter(domain == "faces" & facet == "taste_typicality", sample == "Germine et al. 2015" ) %>% pull(SC_SD)

C_val_sce_eb  = output_C %>% filter(domain == "scenes" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" ) %>% pull(SC_SD)
C_val_fac_eb  = output_C %>% filter(domain == "faces" & facet == "evaluation_bias", sample == "Sutherland et al. 2020" ) %>% pull(SC_SD)

C_val_sce_tt  = output_C %>% filter(domain == "scenes" & facet == "taste_typicality", sample == "Sutherland et al. 2020" ) %>% pull(SC_SD)

#Final Power-analysis:
#Merge (merging simulated and observed power curve)
powCurve_dis_sce_eb_fin = rbind(powCurve_dis_sce_eb %>% mutate(type = "generated"),cbind(powCurve_dis_sce_eb_obs,C = round(C_dis_sce_eb,2)) %>% mutate(type = "observed"))
powCurve_dis_fac_eb_fin = rbind(powCurve_dis_fac_eb %>% mutate(type = "generated"),cbind(powCurve_dis_fac_eb_obs,C = round(C_dis_fac_eb,2)) %>% mutate(type = "observed"))

powCurve_dis_abs_tt_fin = rbind(powCurve_dis_abs_tt %>% mutate(type = "generated"),cbind(powCurve_dis_abs_tt_obs,C = round(C_dis_abs_tt,2)) %>% mutate(type = "observed"))
powCurve_dis_sce_tt_fin = rbind(powCurve_dis_sce_tt %>% mutate(type = "generated"),cbind(powCurve_dis_sce_tt_obs,C = round(C_dis_sce_tt,2)) %>% mutate(type = "observed"))
powCurve_dis_fac_tt_fin = rbind(powCurve_dis_fac_tt %>% mutate(type = "generated"),cbind(powCurve_dis_fac_tt_obs,C = round(C_dis_fac_tt,2)) %>% mutate(type = "observed"))

powCurve_val_sce_eb_fin = rbind(powCurve_val_sce_eb %>% mutate(type = "generated"),cbind(powCurve_val_sce_eb_obs,C = round(C_val_sce_eb,2)) %>% mutate(type = "observed"))
powCurve_val_fac_eb_fin = rbind(powCurve_val_fac_eb %>% mutate(type = "generated"),cbind(powCurve_val_fac_eb_obs,C = round(C_val_fac_eb,2)) %>% mutate(type = "observed"))

powCurve_val_sce_tt_fin = rbind(powCurve_val_sce_tt,cbind(powCurve_val_sce_tt_obs,C = round(C_val_sce_tt,2)))


#plot final power curves
p_pow_dis_sce_eb =   
  ggplot(powCurve_dis_sce_eb_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Evaluation-Bias, Scenes")+
  scale_color_viridis_d()

p_pow_dis_fac_eb =   
  ggplot(powCurve_dis_fac_eb_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Evaluation-Bias, Faces")+
  scale_color_viridis_d()

p_pow_dis_abs_tt =   
  ggplot(powCurve_dis_abs_tt_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Taste-Typicality, Abstracts")+
  scale_color_viridis_d()

p_pow_dis_sce_tt =   
  ggplot(powCurve_dis_sce_tt_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Taste-Typicality, Scenes")+
  scale_color_viridis_d()

p_pow_dis_fac_tt =   
  ggplot(powCurve_dis_fac_tt_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Taste-Typicality, Faces")+
  scale_color_viridis_d()

p_pow_val_sce_eb =   
  ggplot(powCurve_val_sce_eb_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Evaluation-Bias, Scenes")+
  scale_color_viridis_d()

p_pow_val_fac_eb =   
  ggplot(powCurve_val_fac_eb_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C), alpha = type))+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Evaluation-Bias, Faces")+
  scale_color_viridis_d()

p_pow_val_sce_tt =   
  ggplot(powCurve_val_sce_tt_fin %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C)),alpha = .1)+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80,linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis, linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C",
       subtitle = "Taste-Typicality, Scenes")+
  scale_color_viridis_d()


#SAVE####
# pdf(sprintf("%s/05_images/image/supplementary_REV_figure8_ACE_power_disc.pdf",wdNOA),
#     width = 10,
#     height =6 )
# ((plot_spacer()|p_pow_dis_sce_eb|p_pow_dis_fac_eb) /
#     (p_pow_dis_abs_tt|p_pow_dis_sce_tt|p_pow_dis_fac_tt)) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
# dev.off()

png(sprintf("%s/%s/supplementary/04_SF10_ACE_power_discovery.png",wdOA,wdOA_ImageOutput),
    width = 10*325,
    height =6*325,
    unit = "px",
    res = 300,
    bg = "white")
((plot_spacer()|p_pow_dis_sce_eb|p_pow_dis_fac_eb) /
    (p_pow_dis_abs_tt|p_pow_dis_sce_tt|p_pow_dis_fac_tt)) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
dev.off()

png(sprintf("%s/%s/supplementary/04_SF11_ACE_power_validation.png",wdOA,wdOA_ImageOutput),
    width = 5*325,
    height =5*325,
    unit = "px",
    res = 300,
    bg = "white")
((p_pow_val_sce_eb|p_pow_val_fac_eb) /
    (p_pow_val_sce_tt|plot_spacer())) + plot_layout(guides = "collect")+ plot_annotation(tag_levels = "a") & theme(legend.position = "none")
dev.off()


p_pow_area =   
  ggplot(powCurve_dis_fac_eb %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C)),size = 1.5, alpha = .5)+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80, linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis,  linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C")+
  scale_color_viridis_d()


p_pow_area =   
  ggplot(powCurve_dis %>% filter(component == "C"))+
  geom_line(aes(family_n, power, color = as.factor(C)),size = 1.5, alpha = .5)+ 
  geom_point(data=xC_dis,aes(family_n, power))+
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80, linetype = "dashed")+
  geom_vline(xintercept = x_sample_dis,  linetype = "dashed", color = "gray")+
  labs(y = "Power",
       x = "Family-wise sample",
       color = "C")+
  scale_color_viridis_d()
                                      

                                 
p_pow_volm =  powCurve_volm %>% 
  ggplot(aes(family_n, power, color = component,))+
  geom_line(linewidth = 2, alpha = .5)+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80, linetype = "dashed")+
  geom_vline(xintercept = x_sample,  linetype = "dashed", color = "gray")+
  geom_vline(xintercept = xA_volm, linetype = "dashed")+
  geom_vline(xintercept = xC_volm, linetype = "dashed")+
  labs(y = "",
       x = "Family-wise sample",
       subtitle = "Volume")+
  scale_color_manual(values = c("#026A37",
                                         "#538BB5"))
                                         
p_pow_thic =  powCurve_thic %>% 
  ggplot(aes(family_n, power, color = component,))+
  geom_line(linewidth = 2, alpha = .5)+ 
  theme_classic(base_size = 12)+ 
  ylim(0,1)+
  geom_hline(yintercept = .80, linetype = "dashed")+
  geom_vline(xintercept = x_sample, linetype = "dashed", color = "gray")+
  geom_vline(xintercept = xA_thic, linetype = "dashed")+
  geom_vline(xintercept = xC_thic, linetype = "dashed")+
  labs(y = "",
       x = "Family-wise sample",
       subtitle = "Thickness")+
  scale_color_manual(values = c("#026A37",
                                         "#538BB5"))
                                         

(p_pow_area|p_pow_volm|p_pow_thic) +plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

pdf(sprintf("%s/05_images/image/HG_ACE_power.pdf",wdNOA),
    width = 10,
    height =3 )
(p_pow_area|p_pow_volm|p_pow_thic) +plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "none")
dev.off()
