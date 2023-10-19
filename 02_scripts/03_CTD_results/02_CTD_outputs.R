#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 18-09-2023
#Description:
#collate CTD results in two comprhensive outputs: 
#1 Saturated Model compariosn; 
#2 ACDE models
#clean working enviroment 
rm(list = ls())
library(tidyverse)
library(OpenMx)
library(psych)
library(readr)

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#SATURATED####
#load the model comparison test statistics
mod_com_faces_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/07_1_SAT_eb_faces_modelComparison.csv", wdOA,wdOA_output))
mod_com_scences_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/07_2_SAT_eb_scenes_modelComparison.csv", wdOA,wdOA_output))
mod_com_abstracts_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/07_3_SAT_eb_abstracts_modelComparison.csv", wdOA,wdOA_output))

mod_com_faces_tt_dis = read_csv(sprintf("%s/%s/01_Germine_2015/08_1_SAT_tt_faces_modelComparison.csv", wdOA,wdOA_output))
mod_com_scences_tt_dis = read_csv(sprintf("%s/%s/01_Germine_2015/08_2_SAT_tt_scenes_modelComparison.csv", wdOA,wdOA_output))
mod_com_abstracts_tt_dis = read_csv(sprintf("%s/%s/01_Germine_2015/08_3_SAT_tt_abstracts_modelComparison.csv", wdOA,wdOA_output))

mod_com_faces_eb_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/07_1_SAT_eb_faces_modelComparison_val.csv", wdOA,wdOA_output))
mod_com_scences_eb_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/07_2_SAT_eb_scenes_modelComparison_val.csv", wdOA,wdOA_output))

mod_com_faces_tt_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_modelComparison.csv", wdOA,wdOA_output))
mod_com_scences_tt_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/08_2_SAT_tt_scenes_modelComparison.csv", wdOA,wdOA_output))

#bind to a unique df
mod_com = rbind(mod_com_faces_eb_dis,mod_com_scences_eb_dis,mod_com_abstracts_eb_dis,
      mod_com_faces_tt_dis,mod_com_scences_tt_dis,mod_com_abstracts_tt_dis,
      mod_com_faces_eb_val,mod_com_scences_eb_val,
      mod_com_faces_tt_val,mod_com_scences_tt_val)

#tidy for pubblication
mod_com_final = mod_com %>% mutate(base = "Saturated",
                   comparison = recode(comparison, oneCOVAgeca = "Cov=Age", 
                          oneCOVSexca  = "Cov=Sex",
                          oneEMOca = "Mean:Birth_Order",
                          oneEMVOca = "Var:Birth_Order",
                          oneEMVZca = "Mean&Var:Zyg"),
                   comparison = ifelse(is.na(comparison), "Saturated",comparison),
                   minus2LL =round(minus2LL,3) ,
                   AIC = round(AIC,3),
                   diffLL = round(diffLL,3),
                   p = round(p,4),
                   fit = round(fit,3)) %>% 
  select(-c(fitUnits,diffFit,chisq,SBchisq))

#save for supplementary
write_csv(mod_com_final,sprintf("%s/%s/03_CTD_results/02_SF1_SAT_results.csv",wdOA,wdOA_output))

#AC(D)E####
#Do it for CTD ACDE models
ACE_faces_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/09_1_ACE_eb_faces_vc_modelComparison.csv", wdOA,wdOA_output)) %>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ACE_scenes_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/09_2_ACE_eb_scenes_vc_modelComparison.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ADE_abstracts_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/09_3_ADE_eb_abstracts_vc_modelComparsion.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VD",VE =  "VE" ,SA=  "SA",SC_SD = "SD" ,SE= "SE" )

ACE_faces_tt_dis = read_csv(sprintf("%s/%s/01_Germine_2015/10_1_ACE_tt_faces_vc_modelComparsion.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ACE_scenes_tt_dis = read_csv(sprintf("%s/%s/01_Germine_2015/10_2_ACE_tt_scenes_vc_modelComparison.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ACE_abstracts_eb_dis = read_csv(sprintf("%s/%s/01_Germine_2015/10_3_ACE_tt_abstracts_vc_modelComparsion.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )

ACE_faces_eb_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/11_1_ADE_eb_faces_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ACE_scenes_eb_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/11_2_ADE_eb_scenes_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )

ADE_faces_tt_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/12_1_ADE_tt_faces_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VD",VE =  "VE" ,SA=  "SA",SC_SD = "SD" ,SE= "SE" )
ACE_scenes_tt_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/12_2_ACE_tt_scenes_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )

#bind to a unique df
mod_com_ACDE = rbind(ACE_faces_eb_dis,ACE_scenes_eb_dis,ADE_abstracts_eb_dis,
                     ACE_faces_tt_dis,ACE_scenes_tt_dis,ACE_abstracts_eb_dis,
                     ACE_faces_eb_val,ACE_scenes_eb_val,
                     ADE_faces_tt_val,ACE_scenes_tt_val)

mod_com_ACDE_final = mod_com_ACDE %>% mutate(base = recode(base, oneSATca  = "Saturated", 
                                           oneACEvca  = "ACE", 
                                          oneADEvca  = "ADE"),
                             comparison = recode(comparison, oneACEvca  = "ACE", 
                                            oneADEvca  = "ADE", 
                                            oneAEvca  = "AE",
                                            oneCEvca = "CE",
                                            oneEvca = "E"),
                   comparison = ifelse(is.na(comparison), "Saturated",comparison),
                   minus2LL =round(minus2LL,3) ,
                   AIC = round(AIC,3),
                   diffLL = round(diffLL,3),
                   p = round(p,4),
                   fit = round(fit,3)) %>% 
  select(-c(fitUnits,diffFit,chisq,SBchisq))

#save for supplementary
write_csv(mod_com_ACDE_final,sprintf("%s/%s/03_CTD_results/02_SF2_ACE_results.csv",wdOA,wdOA_output))

#Do it for CTD ACDE confidence intervals
ACE_faces_eb_dis_est = read_csv(sprintf("%s/%s/01_Germine_2015/09_1_ACE_eb_faces_vc_bestModel.csv",wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Germine et al. 2015", domain = "faces", facet= "evaluation_bias")
ACE_scenes_eb_dis_est = read_csv(sprintf("%s/%s/01_Germine_2015/09_2_ACE_eb_scenes_vc_bestModel.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Germine et al. 2015", domain = "scenes", facet= "evaluation_bias")
ADE_abstracts_eb_dis_est = read_csv(sprintf("%s/%s/01_Germine_2015/09_3_ADE_eb_abstracts_vc_bestModel.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VD","VE","SA","SD","SE"), sample ="Germine et al. 2015", domain = "abstracts", facet= "evaluation_bias")

ACE_faces_tt_dis_est = read_csv(sprintf("%s/%s/01_Germine_2015/10_1_ACE_tt_faces_vc_bestModel.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Germine et al. 2015", domain = "faces", facet= "taste_typicality")
ACE_scenes_tt_dis_est = read_csv(sprintf("%s/%s/01_Germine_2015/10_2_ACE_tt_scenes_vc_bestModel.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Germine et al. 2015", domain = "scenes", facet= "taste_typicality")
ACE_abstracts_eb_dis_est = read_csv(sprintf("%s/%s/01_Germine_2015/10_3_ACE_tt_abstracts_vc_bestModel.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"),sample ="Germine et al. 2015", domain = "abstracts", facet= "taste_typicality")

ACE_faces_eb_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/11_1_ADE_eb_faces_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Sutherland et al. 2023", domain = "faces", facet= "evaluation_bias")
ACE_scenes_eb_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/11_2_ADE_eb_scenes_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Sutherland et al. 2023", domain = "abstracts", facet= "evaluation_bias")

ADE_faces_tt_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/12_1_ADE_tt_faces_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VD","VE","SA","SD","SE"), sample ="Sutherland et al. 2023", domain = "faces", facet= "taste_typicality")
ACE_scenes_tt_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/12_2_ACE_tt_scenes_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Sutherland et al. 2023", domain = "abstracts", facet= "taste_typicality")

#bind to a unique df
est_ci_ACDE = rbind(ACE_faces_eb_dis_est,ACE_scenes_eb_dis_est,ADE_abstracts_eb_dis_est,
                     ACE_faces_tt_dis_est,ACE_scenes_tt_dis_est,ACE_abstracts_eb_dis_est,
                     ACE_faces_eb_val_est,ACE_scenes_eb_val_est,
                     ADE_faces_tt_val_est,ACE_scenes_tt_val_est)

#save for supplementary
write_csv(est_ci_ACDE,sprintf("%s/%s/03_CTD_results/02_SF3_ACE_ci.csv",wdOA,wdOA_output))


#Table S2####
#estimates and confidence intervals
ACE_faces_ebr_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/13_1_ADE_ebres_faces_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Sutherland et al. 2023", domain = "faces", facet= "evaluation_bias res.")
ACE_scenes_ebr_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/13_2_ADE_ebres_scenes_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Sutherland et al. 2023", domain = "abstracts", facet= "evaluation_bias res.")

ADE_faces_ttr_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/14_1_ADE_ttres_faces_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VD","VE","SA","SD","SE"), sample ="Sutherland et al. 2023", domain = "faces", facet= "taste_typicality res.")
ADE_scenes_ttr_val_est = read_csv(sprintf("%s/%s/02_Sutherland_2020/14_2_ADE_ttres_scenes_vc_bestModel_val.csv", wdOA,wdOA_output))%>% mutate(component = c("VA","VC","VE","SA","SC","SE"), sample ="Sutherland et al. 2023", domain = "abstracts", facet= "taste_typicality res.")

ACE_res =
  rbind(
    #Taste typicality residuals
    ADE_scenes_ttr_val_est,
        ADE_faces_ttr_val_est,
    #Evaluation bias residuals
        ACE_scenes_ebr_val_est,
        ACE_faces_ebr_val_est
  )

View(ACE_res)

#estimates and model comparisons
ACE_faces_ebr_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/13_1_ADE_ebres_faces_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ACE_scenes_ebr_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/13_2_ADE_ebres_scenes_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VC",VE =  "VE" ,SA=  "SA",SC_SD = "SC" ,SE= "SE" )
ADE_faces_ttr_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/14_1_ADE_ttres_faces_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VD",VE =  "VE" ,SA=  "SA",SC_SD = "SD" ,SE= "SE" )
ACE_scenes_ttr_val = read_csv(sprintf("%s/%s/02_Sutherland_2020/14_2_ADE_ttres_scenes_vc_modelComparison_val.csv", wdOA,wdOA_output))%>% rename(VA = "VA", VC_VD = "VD",VE =  "VE" ,SA=  "SA",SC_SD = "SD" ,SE= "SE" )

ACE_res_modelComp = 
  rbind(
ACE_scenes_ttr_val %>% filter(comparison == "oneAEvca") %>% select(minus2LL,AIC,df,diffLL,sample,domain,facet),
ADE_faces_ttr_val %>% filter(comparison == "oneAEvca") %>% select(minus2LL,AIC,df,diffLL,sample,domain,facet),
ACE_scenes_ebr_val %>% filter(comparison == "oneAEvca") %>% select(minus2LL,AIC,df,diffLL,sample,domain,facet),
ACE_faces_ebr_val %>% filter(comparison == "oneAEvca") %>% select(minus2LL,AIC,df,diffLL,sample,domain,facet)
)

View(ACE_res_modelComp)
