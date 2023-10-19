#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 19-09-2023
#
#
#Description: repeat procedure in highlighted in 01_Germine_2015/08_1_Fac_mm2_UnivariateSAT.R
#clean working environment 
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
BioMetric_Pheno  = read_csv(sprintf("%s/%s/02_Sutherland_2020/05_twin_tt.csv", wdOA,wdOA_output))

#load functions:
source(sprintf("%s/%s/functions/miFunctions.R", wdOA,wdOA_scripts))

Fac_BioMetric_Pheno = BioMetric_Pheno%>%
  dplyr::select(Sub, AgeAtTestingCalculated, Sex, FA_A,FA_B,TwinType)%>%
  rename(Age = "AgeAtTestingCalculated", 
         mm2_z_FA_1 = "FA_A",
         mm2_z_FA_2 = "FA_B")%>%
  as.data.frame()



# REPARE DATA ----------------------------------------------------------------------------------------------------------------------
# Select Variables for Analysis
vars <- 'mm2_z_FA' # list of variables names
nv <- 1 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")
# Select Covariates for Analysis
covVars <- c('Age', 'Sex')
# Select Data for Analysis
mzData <- subset(Fac_BioMetric_Pheno, TwinType==1, c(selVars, covVars))
dzData <- subset(Fac_BioMetric_Pheno, TwinType==3, c(selVars, covVars))
# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Set Starting Values
svBeAge <- -0.01 # start value for regressions
svBeSex <- 0.01 # start value for regressions
svMe <- .9 # start value for means
svVa <- .2 # start value for variance
# lbVa <- .02 # lower bound for variance
# PREPARE MODEL----------------------------------------------------------------------------------------------------------------------
# Create Matrices for Covariates and linear Regression Coefficients
defLage <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.Age"), name="defLage" )
defLsex <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.Sex"), name="defLsex" )
pathBlage <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBeAge, label="b11", name="blage" )
pathBlsex <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBeSex, label="b12", name="blsex" )

# Create Algebra for expected Mean Matrices
meanMZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZ1","mMZ2"), name="meanMZ" )
meanDZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZ1","mDZ2"), name="meanDZ" )
expMeanMZ <- mxAlgebra( expression= meanMZ + cbind(defLage%*%blage,defLage%*%blage) + cbind(defLsex%*%blsex,defLsex%*%blsex), name="expMeanMZ" )
expMeanDZ <- mxAlgebra( expression= meanDZ + cbind(defLage%*%blage,defLage%*%blage) + cbind(defLsex%*%blsex,defLsex%*%blsex), name="expMeanDZ" )
# Create Algebra for expected Variance/Covariance Matrices
covMZ <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=c("vMZ1","cMZ21","vMZ2"),#lbound=valDiag(lbVa,ntv)
                   name="covMZ" )
covDZ <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=c("vDZ1","cDZ21","vDZ2"),#lbound=valDiag(lbVa,ntv)
                   name="covDZ" )
# Calculate correlations from expected Covariance Matrices
corMZ      <- mxAlgebra( expression=cov2cor(covMZ), name ="corMZ")
corDZ     <- mxAlgebra( expression=cov2cor(covDZ), name ="corDZ")

# Create Data Objects for Multiple Groups
dataMZ <- mxData( observed=mzData, type="raw" )
dataDZ <- mxData( observed=dzData, type="raw" )
# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal( covariance="covMZ", means="expMeanMZ", dimnames=selVars )
expDZ <- mxExpectationNormal( covariance="covDZ", means="expMeanDZ", dimnames=selVars )
funML <- mxFitFunctionML()
# Create Model Objects for Multiple Groups
pars <- list( pathBlage, pathBlsex )
defs <- list( defLage, defLsex )
modelMZ <- mxModel( pars, defs, meanMZ, expMeanMZ, covMZ, corMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ <- mxModel( pars, defs, meanDZ, expMeanDZ, covDZ, corDZ, dataDZ, expDZ, funML, name="DZ" )
multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
# Create Confidence Interval Objects
ciCov <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciCor <- mxCI( c('MZ.corMZ[2,1]','DZ.corDZ[2,1]') )
ciMean <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )
# Build Saturated Model with Confidence Intervals
modelSAT <- mxModel( "oneSATca", pars, modelMZ, modelDZ, multi, ciCov, ciCor,ciMean )
# RUN MODEL----------------------------------------------------------------------------------------------------------------------
####_Saturated Model####
# Run Saturated Model
fitSAT <- mxRun( modelSAT, intervals=F )
sumSAT <- summary( fitSAT )
# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitSAT)
fitEsts(fitSAT)
mxGetExpected( fitSAT, c("means","covariance") )

# RUN SUBMODELS----------------------------------------------------------------------------------------------------------------------
####_Covariate####
# Test Significance of Covariate: age
modelCOVage <- mxModel( fitSAT, name="oneCOVAgeca" )
modelCOVage <- omxSetParameters( modelCOVage, label="b11", free=FALSE, values=0 )
fitCOVage <- mxRun( modelCOVage )
# Test Significance of Covariate: age
modelCOVsex <- mxModel( fitSAT, name="oneCOVSexca" )
modelCOVsex <- omxSetParameters( modelCOVsex, label="b12", free=FALSE, values=0 )
fitCOVsex <- mxRun( modelCOVsex )
####_Birth Order####
# Constrain expected Means to be equal across Twin Order
modelEMO <- mxModel( fitSAT, name="oneEMOca" )
modelEMO <- omxSetParameters( modelEMO, label=c("mMZ1","mMZ2"), free=TRUE, values=svMe, newlabels='mMZ' )
modelEMO <- omxSetParameters( modelEMO, label=c("mDZ1","mDZ2"), free=TRUE, values=svMe, newlabels='mDZ' )
fitEMO <- mxRun( modelEMO, intervals=F )
fitGofs(fitEMO); fitEsts(fitEMO)
# Constrain expected Means and Variances to be equal across Twin Order
modelEMVO <- mxModel( fitEMO, name="oneEMVOca" )
modelEMVO <- omxSetParameters( modelEMVO, label=c("vMZ1","vMZ2"), free=TRUE, values=svVa, newlabels='vMZ' )
modelEMVO <- omxSetParameters( modelEMVO, label=c("vDZ1","vDZ2"), free=TRUE, values=svVa, newlabels='vDZ' )
fitEMVO <- mxRun( modelEMVO, intervals=F )
fitGofs(fitEMVO); fitEsts(fitEMVO)
####_Zygosity####
# Constrain expected Means and Variances to be equal across Twin Order and Zygosity
modelEMVZ <- mxModel( fitEMVO, name="oneEMVZca" )
modelEMVZ <- omxSetParameters( modelEMVZ, label=c("mMZ","mDZ"), free=TRUE, values=svMe, newlabels='mZ' )
modelEMVZ <- omxSetParameters( modelEMVZ, label=c("vMZ","vDZ"), free=TRUE, values=svVa, newlabels='vZ' )
fitEMVZ <- mxRun( modelEMVZ, intervals=T )
fitGofs(fitEMVZ); fitEsts(fitEMVZ)
####COMPARISON####
# Print Comparative Fit Statistics
modelComparison = mxCompare( fitSAT, subs <- list(fitCOVage, fitCOVsex, fitEMO, fitEMVO, fitEMVZ) )%>%
  as.data.frame()%>%
  mutate(sample = "Sutherland et al. 2020",
         domain = "faces",
         facet = "taste_typicality")

# PARAMETERS----------------------------------------------------------------------------------------------------------------------
####_mean and standard deviation####
fitSATdescriptives = fitEstCis(fitEMVZ)%>%
  as.data.frame()%>%
  rownames_to_column()%>%
  rename(parameter = "rowname")

fitSATdescriptives[fitSATdescriptives$parameter == "mZ",]$estimate#mean
####_expected model####

fitSATdescriptives$expected_full = ifelse( fitSATdescriptives[fitSATdescriptives$parameter == "DZ.corDZ[2,1]", ]$estimate >.5*(fitSATdescriptives[fitSATdescriptives$parameter == "MZ.corMZ[2,1]", ]$estimate), "ACE", "ADE" )
fitSATdescriptives$expected_reduced =  ifelse( fitSATdescriptives[fitSATdescriptives$parameter == "DZ.corDZ[2,1]", ]$estimate < ( fitSATdescriptives[fitSATdescriptives$parameter == "MZ.corMZ[2,1]", ]$estimate ), "AE", "CE" )

# SAVE----------------------------------------------------------------------------------------------------------------------
save(fitSAT,fitCOVsex,fitCOVsex,fitEMO,fitEMO,fitEMVO,fitEMVZ, file = sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_Sutherland2020.RData",wdOA,wdOA_output))
save(modelComparison,file = sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_modelComparison.RData",wdOA,wdOA_output))
write_csv(fitSATdescriptives,sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_parameters.csv",wdOA,wdOA_output))
write_csv(modelComparison,sprintf("%s/%s/02_Sutherland_2020/08_1_SAT_tt_faces_modelComparison.csv",wdOA,wdOA_output))
