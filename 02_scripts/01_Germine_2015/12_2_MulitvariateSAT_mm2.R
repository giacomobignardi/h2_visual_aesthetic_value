#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
# Twin Multivariate (trivariate) Saturated Structural Equation Model: 2 groups (MZss and DZss)
# Matrix style model - Raw data - Continuous data
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

#load dataFrames:
BioMetric_Pheno  = read_csv(sprintf("%s/%s/01_Germine_2015/05_twin_tt_multivariate.csv", wdOA,wdOA_output))

#load functions:
source(sprintf("%s/%s/functions/miFunctions.R", wdOA,wdOA_scripts))

muvTraits = c("AO_1","AO_2","FA_TOT_1","FA_TOT_2", "SC_1", "SC_2")

BioMetric_multivariate = BioMetric_Pheno%>%
  select(c(FamId,Zygosity),muvTraits)
describe(BioMetric_multivariate[,c("AO_1","AO_2","FA_TOT_1","FA_TOT_2", "SC_1", "SC_2")], skew=F)

# Select Variables for Analysis
vars      <- c('SC','FA_TOT','AO') 
nv        <- 3       # number of variables
ntv       <- nv*2    # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")

#Select Data for Analysis
mzData    <- subset(BioMetric_multivariate, Zygosity==1, selVars)
dzData    <- subset(BioMetric_multivariate, Zygosity==3, selVars)

# Generate Descriptive Statistics
round(colMeans(mzData,na.rm=TRUE),4)
round(colMeans(dzData,na.rm=TRUE),4)
round(cov(mzData,use="complete"),4)
round(cov(dzData,use="complete"),4)
round(cor(mzData,use="complete"),4)
round(cor(dzData,use="complete"),4)

# Set Starting Values 
svMe      <- c(0,0,0)        # start value for means
svVa      <- 1                     # start value for variances
lbVa      <- .01                 # start value for lower bounds

# Create Labels
labMeMZ   <- paste("meanMZ",selVars,sep="_")
labMeDZ   <- paste("meanDZ",selVars,sep="_")
labMeZ    <- paste("meanZ",selVars,sep="_")
labCvMZ   <- labLower("covMZ",ntv)
labCvDZ   <- labLower("covDZ",ntv)
labCvZ    <- labLower("covZ",ntv)
labVaMZ   <- labDiag("covMZ",ntv)
labVaDZ   <- labDiag("covDZ",ntv)
labVaZ    <- labDiag("covZ",ntv)

# ------------------------------------------------------------------------------
# PREPARE MODEL
set.seed(42)
# Saturated Model
# Create Algebra for expected Mean Matrices
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeMZ, name="meanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZ, name="meanDZ" )

# Create Algebra for expected Variance/Covariance Matrices
covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=labCvMZ, name="covMZ" )
covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=labCvDZ, name="covDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
modelMZ   <- mxModel( "MZ", meanMZ, covMZ, dataMZ, expMZ, funML )
modelDZ   <- mxModel( "DZ", meanDZ, covDZ, dataDZ, expDZ, funML )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Confidence Interval Objects
ciCov     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )

# Build Saturated Model with Confidence Intervals
modelSAT  <- mxModel( "mulSATc", modelMZ, modelDZ, multi, ciCov, ciMean )

# ------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
fitSAT = mxRun( modelSAT, intervals=T )
# fit       <- mxTryHard( fitSAT )
sumSAT  = summary( fitSAT )
fitGofs(fitSAT)
fitEstCis(fitSAT)

#Parameters------------------------------------------------------------------------------     
Param_fit_means = mxGetExpected(fitSAT, c("means"))
Param_fit_cov = mxGetExpected(fitSAT, c("covariance"))

MZmeans = Param_fit_means$MZ
MZcovCor = as.data.frame(cov2cor(Param_fit_cov$MZ))
DZmeans = Param_fit_means$DZ 
DZcovCor = as.data.frame(cov2cor(Param_fit_cov$DZ))

####SAVE####
save(fitSAT, file = sprintf("%s/%s/01_Germine_2015/12_2_SAT_tt_multivariate_Germine2015.RData",wdOA,wdOA_output))
write_csv(MZcovCor%>%rownames_to_column(),sprintf("%s/%s/01_Germine_2015/12_2_SAT_tt_MZcorMatrix_multivariate.csv",wdOA,wdOA_output))
write_csv(DZcovCor%>%rownames_to_column(),sprintf("%s/%s/01_Germine_2015/12_2_SAT_tt_DZcorMatrix_multivariate.csv",wdOA,wdOA_output))


