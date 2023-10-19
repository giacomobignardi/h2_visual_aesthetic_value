#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description: repeat VCA (as in 02_mlm_twin1) for twin 2
#Program: 02_supplementary_mlm_twin" ------------------------------------------------------------------------------------------------------------------------------

#clear workspace
rm(list = ls())

#load packages
library(tidyverse)
library(lme4)
library(pbkrtest)

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#load dataFrames:
Fac_MLM_TwinLong  = read_csv(sprintf("%s/%s/01_Germine_2015/01_MLM_Fac_Twin2_Germine2015.csv", wdOA,wdOA_output))
Sce_MLM_TwinLong  = read_csv(sprintf("%s/%s/01_Germine_2015/01_MLM_Sce_Twin2_Germine2015.csv", wdOA,wdOA_output))
Obj_MLM_TwinLong  = read_csv(sprintf("%s/%s/01_Germine_2015/01_MLM_Obj_Twin2_Germine2015.csv", wdOA,wdOA_output))
#load functions:
source(sprintf("%s/%s/functions/VCA.R", wdOA,wdOA_scripts))

#parameters:
#simulations (no need to compute CI here)
NSim = 1
#N simulation for the semi-parametric bootstraping (NB if a model takes time=t to fit, n simulations will take t=n*t)
set.seed(42)

#rename first 4 columns to fit entry requirements for VCA function
names(Fac_MLM_TwinLong) = c("Obj","Rating","Sub","Block")
names(Sce_MLM_TwinLong) = c("Obj","Rating","Sub","Block")
names(Obj_MLM_TwinLong) = c("Obj","Rating","Sub","Block")

####MLM####
old <- Sys.time() # get start time
#fit multilevel models (based on Martinez et al 2020)
Opt = lmerControl(optimizer = "bobyqa", calc.derivs = F)
Fac_MLM <- lmer(Rating ~ 1 + ((1 | Sub) + (1 | Obj) + (1|Sub:Obj) + (1|Block) + (1|Block:Sub) + (1|Block:Obj)), data=Fac_MLM_TwinLong, control = Opt)
Sce_MLM <- lmer(Rating ~ 1 + ((1 | Sub) + (1 | Obj) + (1|Sub:Obj) + (1|Block) + (1|Block:Sub) + (1|Block:Obj)), data=Sce_MLM_TwinLong, control = Opt)
Obj_MLM <- lmer(Rating ~ 1 + ((1 | Sub) + (1 | Obj) + (1|Sub:Obj) + (1|Block) + (1|Block:Sub) + (1|Block:Obj)), data=Obj_MLM_TwinLong, control = Opt)

#semi parametric bootstraping to estimate CIs
options(nwarnings = 500) #set warning to see how many model did not converged during the simulations
Fac_bootVCA = bootMer(Fac_MLM, VCA, nsim=NSim, .progress = "txt")
Fac_boot_Warn = warnings()
Sce_bootVCA = bootMer(Sce_MLM, VCA, nsim=NSim, .progress = "txt")
Sce_boot_Warn = warnings()
Obj_bootVCA = bootMer(Obj_MLM, VCA, nsim=NSim, .progress = "txt")
Obj_boot_Warn = warnings()
# print elapsed time
new <- Sys.time() - old # calculate difference
print(old)
print(new)
#etimated time of completion 9168.425m (~6d)

####VCA####
####_faces####
#Create a summuary of the VCA
Fac_VCAsummary = 
  data.frame(
  Domain = rep("FA_TOT",7),
  Component = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock", "Residual"),
  Value = VCA(Fac_MLM),
  CI_low = c(
    as.numeric(quantile(Fac_bootVCA$t[,1], c(0.025, 0.975)))[1],
    as.numeric(quantile(Fac_bootVCA$t[,2], c(0.025, 0.975)))[1],
    as.numeric(quantile(Fac_bootVCA$t[,3], c(0.025, 0.975)))[1],
    as.numeric(quantile(Fac_bootVCA$t[,4], c(0.025, 0.975)))[1],
    as.numeric(quantile(Fac_bootVCA$t[,5], c(0.025, 0.975)))[1],
    as.numeric(quantile(Fac_bootVCA$t[,6], c(0.025, 0.975)))[1],
    as.numeric(quantile(Fac_bootVCA$t[,7], c(0.025, 0.975)))[1]
  ),
  CI_high = c(
    as.numeric(quantile(Fac_bootVCA$t[,1], c(0.025, 0.975)))[2],
    as.numeric(quantile(Fac_bootVCA$t[,2], c(0.025, 0.975)))[2],
    as.numeric(quantile(Fac_bootVCA$t[,3], c(0.025, 0.975)))[2],
    as.numeric(quantile(Fac_bootVCA$t[,4], c(0.025, 0.975)))[2],
    as.numeric(quantile(Fac_bootVCA$t[,5], c(0.025, 0.975)))[2],
    as.numeric(quantile(Fac_bootVCA$t[,6], c(0.025, 0.975)))[2],
    as.numeric(quantile(Fac_bootVCA$t[,7], c(0.025, 0.975)))[2]
  )
  )

####_sceneries####
#Create a summuary of the VCA
Sce_VCAsummary = 
  data.frame(
    Domain = rep("SC",7),
    Component = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock", "Residual"),
    Value = VCA(Sce_MLM),
    CI_low = c(
      as.numeric(quantile(Sce_bootVCA$t[,1], c(0.025, 0.975)))[1],
      as.numeric(quantile(Sce_bootVCA$t[,2], c(0.025, 0.975)))[1],
      as.numeric(quantile(Sce_bootVCA$t[,3], c(0.025, 0.975)))[1],
      as.numeric(quantile(Sce_bootVCA$t[,4], c(0.025, 0.975)))[1],
      as.numeric(quantile(Sce_bootVCA$t[,5], c(0.025, 0.975)))[1],
      as.numeric(quantile(Sce_bootVCA$t[,6], c(0.025, 0.975)))[1],
      as.numeric(quantile(Sce_bootVCA$t[,7], c(0.025, 0.975)))[1]
    ),
    CI_high = c(
      as.numeric(quantile(Sce_bootVCA$t[,1], c(0.025, 0.975)))[2],
      as.numeric(quantile(Sce_bootVCA$t[,2], c(0.025, 0.975)))[2],
      as.numeric(quantile(Sce_bootVCA$t[,3], c(0.025, 0.975)))[2],
      as.numeric(quantile(Sce_bootVCA$t[,4], c(0.025, 0.975)))[2],
      as.numeric(quantile(Sce_bootVCA$t[,5], c(0.025, 0.975)))[2],
      as.numeric(quantile(Sce_bootVCA$t[,6], c(0.025, 0.975)))[2],
      as.numeric(quantile(Sce_bootVCA$t[,7], c(0.025, 0.975)))[2]
    )
  )

####_abstracts####
#Create a summuary of the VCA
Obj_VCAsummary = 
  data.frame(
    Domain = rep("AO",7),
    Component = c("Obj","Sub","Block", "SubXObj", "SubXBlock", "ObjxBlock", "Residual"),
    Value = VCA(Obj_MLM),
    CI_low = c(
      as.numeric(quantile(Obj_bootVCA$t[,1], c(0.025, 0.975)))[1],
      as.numeric(quantile(Obj_bootVCA$t[,2], c(0.025, 0.975)))[1],
      as.numeric(quantile(Obj_bootVCA$t[,3], c(0.025, 0.975)))[1],
      as.numeric(quantile(Obj_bootVCA$t[,4], c(0.025, 0.975)))[1],
      as.numeric(quantile(Obj_bootVCA$t[,5], c(0.025, 0.975)))[1],
      as.numeric(quantile(Obj_bootVCA$t[,6], c(0.025, 0.975)))[1],
      as.numeric(quantile(Obj_bootVCA$t[,7], c(0.025, 0.975)))[1]
    ),
    CI_high = c(
      as.numeric(quantile(Obj_bootVCA$t[,1], c(0.025, 0.975)))[2],
      as.numeric(quantile(Obj_bootVCA$t[,2], c(0.025, 0.975)))[2],
      as.numeric(quantile(Obj_bootVCA$t[,3], c(0.025, 0.975)))[2],
      as.numeric(quantile(Obj_bootVCA$t[,4], c(0.025, 0.975)))[2],
      as.numeric(quantile(Obj_bootVCA$t[,5], c(0.025, 0.975)))[2],
      as.numeric(quantile(Obj_bootVCA$t[,6], c(0.025, 0.975)))[2],
      as.numeric(quantile(Obj_bootVCA$t[,7], c(0.025, 0.975)))[2]
    )
  )


####modified BI####
#Compute Beholder indcies modified (Honnekop, 2006)
#modification is based on the number of component extracted from the MLM

####_faces####
#calculate repearable and error variance
Fac_VCAsummary$Repeatable = (Fac_VCAsummary[Fac_VCAsummary$Component == "Sub",]$Value +
                               Fac_VCAsummary[Fac_VCAsummary$Component == "Obj",]$Value + 
                               Fac_VCAsummary[Fac_VCAsummary$Component == "SubXObj",]$Value + 
                               Fac_VCAsummary[Fac_VCAsummary$Component == "SubXBlock",]$Value + 
                               Fac_VCAsummary[Fac_VCAsummary$Component == "ObjxBlock",]$Value)
Fac_VCAsummary$Error = 1-Fac_VCAsummary$Repeatable

#calculate mBI: assumption Sub component is not error variance
Fac_VCAsummary$mBi = (Fac_VCAsummary[Fac_VCAsummary$Component == "Sub",]$Value +
                        Fac_VCAsummary[Fac_VCAsummary$Component == "SubXObj",]$Value + 
                        Fac_VCAsummary[Fac_VCAsummary$Component == "SubXBlock",]$Value)/
  (Fac_VCAsummary$Repeatable[1])


####_scenes####
#calculate repearable and error variance
Sce_VCAsummary$Repeatable = (Sce_VCAsummary[Sce_VCAsummary$Component == "Sub",]$Value +
                               Sce_VCAsummary[Sce_VCAsummary$Component == "Obj",]$Value + 
                               Sce_VCAsummary[Sce_VCAsummary$Component == "SubXObj",]$Value + 
                               Sce_VCAsummary[Sce_VCAsummary$Component == "SubXBlock",]$Value + 
                               Sce_VCAsummary[Sce_VCAsummary$Component == "ObjxBlock",]$Value)
Sce_VCAsummary$Error = 1-Sce_VCAsummary$Repeatable 

#calculate mBI: assumption Sub component is not error variance
Sce_VCAsummary$mBi = (Sce_VCAsummary[Sce_VCAsummary$Component == "Sub",]$Value +
                        Sce_VCAsummary[Sce_VCAsummary$Component == "SubXObj",]$Value + 
                        Sce_VCAsummary[Sce_VCAsummary$Component == "SubXBlock",]$Value)/
  (Sce_VCAsummary$Repeatable[1])

####_abstracts####
#calculate repearable and error variance
Obj_VCAsummary$Repeatable = (Obj_VCAsummary[Obj_VCAsummary$Component == "Sub",]$Value +
                               Obj_VCAsummary[Obj_VCAsummary$Component == "Obj",]$Value + 
                               Obj_VCAsummary[Obj_VCAsummary$Component == "SubXObj",]$Value + 
                               Obj_VCAsummary[Obj_VCAsummary$Component == "SubXBlock",]$Value + 
                               Obj_VCAsummary[Obj_VCAsummary$Component == "ObjxBlock",]$Value)
Obj_VCAsummary$Error = 1-Obj_VCAsummary$Repeatable

#calculate mBI: assumption Sub component is not error variance
Obj_VCAsummary$mBi = (Obj_VCAsummary[Obj_VCAsummary$Component == "Sub",]$Value +
                        Obj_VCAsummary[Obj_VCAsummary$Component == "SubXObj",]$Value + 
                        Obj_VCAsummary[Obj_VCAsummary$Component == "SubXBlock",]$Value)/
  (Obj_VCAsummary$Repeatable[1])

VCAsummary = rbind(Fac_VCAsummary,Sce_VCAsummary, Obj_VCAsummary)

#####block variance####
#calculate total amount of variance explained by block and its interactions
VCAsummary%>%
  filter(Component == "Block" | Component == "SubXBlock" | Component == "ObjxBlock")%>%
  group_by(Domain)%>%
  summarise(sum(Value))

####SAVE###
MLMsummary =  c(Fac_MLM, Fac_bootVCA, Sce_MLM, Sce_bootVCA, Obj_MLM, Obj_bootVCA)
write_csv(VCAsummary,sprintf("%s/%s/01_Germine_2015/02_sup_VCA_twin2_Germine2015.csv",wdOA,wdOA_output))
save(MLMsummary, file = sprintf("%s/%s/01_Germine_2015/02_sup_MLM_twin2_Germine2015.RData",wdOA,wdOA_output))
