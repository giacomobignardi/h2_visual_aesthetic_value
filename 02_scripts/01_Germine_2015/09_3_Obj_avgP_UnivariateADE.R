#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
# Twin Univariate CTD informed Structural Equation Model to estimate variance components (A,C or D, and E) across 2 groups (MZss and DZss)
# Matrix style model - Raw data - Continuous data
#clean working enviroment 
rm(list = ls())
library(tidyverse)
library(OpenMx)
library(psych)

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#load dataFrames:
BioMetric  = read_csv(sprintf("%s/%s/01_Germine_2015/05_twin_univariate.csv",wdOA,wdOA_output))
#load functions:
source(sprintf("%s/%s/functions/miFunctions.R", wdOA,wdOA_scripts))
set.seed(42)
Obj_BioMetric_Pheno = BioMetric%>%filter(category =="AO")%>%dplyr::select(FamId, Age, Sex, avgP_1,avgP_2,Zygosity)%>%as.data.frame()
# ----------------------------------------------------------------------------------------------------------------------
vars <- 'avgP' # list of variables names
nv <- 1 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")
# Select Covariates for Analysis
covVars <- c('Age','Sex')
# Select Data for Analysis
mzData <- subset(Obj_BioMetric_Pheno, Zygosity==1, c(selVars, covVars))
dzData <- subset(Obj_BioMetric_Pheno, Zygosity==3, c(selVars, covVars))
# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Set Starting Values
svBe <- 0.01 # start value for regressions
svMe <- 3 # start value for means
svPa <- .3 # start value for path coefficient
svPe <- .6 # start value for path coefficient for e
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL
# Create Matrices for Covariates and linear Regression Coefficients
defLage <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.Age"), name="defLage" )
defLsex <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.Sex"), name="defLsex" )
pathBlage <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBe, label="b11", name="blage" )
pathBlsex <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBe, label="b12", name="blsex" )
# Create Algebra for expected Mean Matrices
meanG <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )
expMean <- mxAlgebra( expression= meanG + cbind(defLage%*%blage,defLage%*%blage) + cbind(defLsex%*%blsex,defLsex%*%blsex), name="expMeanG" )
# Create Matrices for Variance Components
covA <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" )
covD <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VD11", name="VD" )
covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VE11", name="VE" )
# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra( expression= VA+VD+VE, name="V" )
covMZ <- mxAlgebra( expression= VA+VD, name="cMZ" )
covDZ <- mxAlgebra( expression= 0.5%x%VA+ 0.25%x%VD, name="cDZ" )
expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
# Create Data Objects for Multiple Groups
dataMZ <- mxData( observed=mzData, type="raw" )
dataDZ <- mxData( observed=dzData, type="raw" )
# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal( covariance="expCovMZ", means="expMeanG", dimnames=selVars )
expDZ <- mxExpectationNormal( covariance="expCovDZ", means="expMeanG", dimnames=selVars )
funML <- mxFitFunctionML()
# Create Model Objects for Multiple Groups
pars <- list( pathBlage, pathBlsex , meanG, covA, covD, covE, covP )
defs <- list( defLage, defLsex )
modelMZ <- mxModel( pars, defs, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
# Create Algebra for Variance Components
rowUS <- rep('US',nv)
colUS <- rep(c('VA','VD','VE','SA','SD','SE'),each=nv)
estUS <- mxAlgebra( expression=cbind(VA,VD,VE,VA/V,VD/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
# Create Confidence Interval Objects
ciADE <- mxCI( "US[1,1:6]" )
# Build Model with Confidence Intervals
modelADE <- mxModel( "oneADEvca", pars, modelMZ, modelDZ, multi, estUS, ciADE )
# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL
# Run ADE Model
fitADE <- mxRun( modelADE, intervals=T )
sumADE <- summary( fitADE )
# Compare with Saturated Model
#if saturated model fitted in same session
#mxCompare( fitSAT, fitACE )
#if saturated model prior to genetic model
#lrtSAT(fitADE,4055.9346,1767)
# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitADE)
fitEstCis(fitADE)
# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS
# Run AE model
modelAE <- mxModel( fitADE, name="oneAEvca" )
modelAE <- omxSetParameters( modelAE, labels="VD11", free=FALSE, values=0 )
fitAE <- mxRun( modelAE, intervals=T )
fitGofs(fitAE); fitEstCis(fitAE)
# Run E model
modelE <- mxModel( fitAE, name="oneEvca" )
modelE <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
fitE <- mxRun( modelE, intervals=T )
fitGofs(fitE); fitEstCis(fitE)
#load dataFrames:
load(sprintf("%s/%s/01_Germine_2015/07_3_SAT_eb_abstracts_Germine2015.RData", wdOA,wdOA_output))
# Print Comparative Fit Statistics
modelComparison = rbind(
  mxCompare( fitSAT, fitADE )%>% as.data.frame(),
  mxCompare( fitADE, nested <- list(fitAE, fitE) ) %>% as.data.frame() 
) %>% 
  mutate(sample = "Germine et al. 2015",
         domain = "abstracts",
         facet = "evaluation_bias")%>% 
  filter(!row_number() %in% c(3))


ADEvc = rbind(rep(NA,6),as.data.frame(round(rbind(fitADE$US$result,fitAE$US$result,fitE$US$result),4)))
#Best model VCAest
AECIests = fitEstCis(fitAE)
fitEstCis(fitADE)
# SAVE----------------------------------------------------------------------------------------------------------------------
save(fitADE,fitAE,fitE,fitE, file = sprintf("%s/%s/01_Germine_2015/09_3_ADE_eb_abstracts_Germine2015.RData",wdOA,wdOA_output))
write_csv(cbind(modelComparison,ADEvc),sprintf("%s/%s/01_Germine_2015/09_3_ADE_eb_abstracts_vc_modelComparsion.csv",wdOA,wdOA_output))
write_csv(as.data.frame(AECIests),sprintf("%s/%s/01_Germine_2015/09_3_ADE_eb_abstracts_vc_bestModel.csv",wdOA,wdOA_output))

