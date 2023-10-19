#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 19-09-2023
#
#
#Description: repeat procedure in highlighted in 01_Germine_2015/XXX.R

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
BioMetric_Pheno  = read_csv(sprintf("%s/%s/02_Sutherland_2020/05_twin_eb.csv", wdOA,wdOA_output))

#load functions:
source(sprintf("%s/%s/functions/miFunctions.R", wdOA,wdOA_scripts))

Fac_BioMetric_Pheno = BioMetric_Pheno%>%
  dplyr::select(Sub, AgeAtTestingCalculated, Sex, FA_A, FA_B,TwinType)%>%
  rename(Age = "AgeAtTestingCalculated",
         avgP_FA_1 = "FA_A",
         avgP_FA_2 = "FA_B")%>%
  as.data.frame()

# ----------------------------------------------------------------------------------------------------------------------
vars <- 'avgP_FA' # list of variables names
nv <- 1 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")
# Select Covariates for Analysis
covVars <- c('Age','Sex')
# Select Data for Analysis
mzData <- subset(Fac_BioMetric_Pheno, TwinType==1, c(selVars, covVars))
dzData <- subset(Fac_BioMetric_Pheno, TwinType==3, c(selVars, covVars))
# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Set Starting Values
svBe <- 0.0001 # start value for regressions
svMe <- 0 # start value for means
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
covC <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VC11", name="VC" )
covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VE11", name="VE" )
# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
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
pars <- list( pathBlage, pathBlsex, meanG, covA, covC, covE, covP )
defs <- list( defLage, defLsex  )
modelMZ <- mxModel( pars, defs, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
# Create Algebra for Variance Components
rowUS <- rep('US',nv)
colUS <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
estUS <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
# Create Confidence Interval Objects
ciACE <- mxCI( "US[1,1:6]" )
# Build Model with Confidence Intervals
modelACE <- mxModel( "oneACEvca", pars, modelMZ, modelDZ, multi, estUS, ciACE )
# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL
# Run ACE Model
fitACE <- mxRun( modelACE, intervals=T )
sumACE <- summary( fitACE )
# Compare with Saturated Model
#if saturated model fitted in same session
#mxCompare( fitSAT, fitACE )
#if saturated model prior to genetic model
#lrtSAT(fitACE,4055.9346,1767)
# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACE)
fitEstCis(fitACE)
# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS
# Run AE model
modelAE <- mxModel( fitACE, name="oneAEvca" )
modelAE <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
fitAE <- mxRun( modelAE, intervals=T )
fitGofs(fitAE); fitEstCis(fitAE)
# Run CE model
modelCE <- mxModel( fitACE, name="oneCEvca" )
modelCE <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
modelCE <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
fitCE <- mxRun( modelCE, intervals=T )
fitGofs(fitCE); fitEstCis(fitCE)
# Run E model
modelE <- mxModel( fitAE, name="oneEvca" )
modelE <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
fitE <- mxRun( modelE, intervals=T )
fitGofs(fitE); fitEstCis(fitE)
# Print Comparative Fit Statistics
load(sprintf("%s/%s/02_Sutherland_2020/07_1_SAT_eb_faces_Sutherland2020.Rdata", wdOA,wdOA_output))
#compare models
modelComparison = rbind(
  mxCompare( fitSAT, fitACE )%>% as.data.frame(),
  mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) ) %>% as.data.frame() 
) %>% 
  mutate(sample = "Sutherland et al. 2020",
         domain = "faces",
         facet = "evaluation_bias") %>% 
  filter(!row_number() %in% c(3))


ACEvc = rbind(rep(NA,6),as.data.frame(round(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4)))

#Best model VCAest
AECIests = fitEstCis(fitAE)
fitEstCis(fitACE)
# SAVE----------------------------------------------------------------------------------------------------------------------
save(fitACE,fitAE,fitCE,fitE, file = sprintf("%s/%s/02_Sutherland_2020/11_1_ADE_eb_faces_Sutherland2020.Rdata",wdOA,wdOA_output))
write_csv(cbind(modelComparison,ACEvc),sprintf("%s/%s/02_Sutherland_2020/11_1_ADE_eb_faces_vc_modelComparison_val.csv",wdOA,wdOA_output))
write_csv(as.data.frame(AECIests),sprintf("%s/%s/02_Sutherland_2020/11_1_ADE_eb_faces_vc_bestModel_val.csv",wdOA,wdOA_output))


