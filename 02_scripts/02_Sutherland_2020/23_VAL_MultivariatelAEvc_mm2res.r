#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#
#
#
#Description:
# Twin Univariate Saturated model to estimate means and (co)variances across multiple groups
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

BioMetric_Pheno  = read_csv(sprintf("%s/%s/02_Sutherland_2020/05_twin_ttres_multivariate.csv", wdOA,wdOA_output))
#load fit sat to compare with saturated model
load(sprintf("%s/%s/02_Sutherland_2020/19_SAT_ttres_multivariate_Sutherland2020.RData", wdOA,wdOA_output))

#load functions:
source(sprintf("%s/%s/functions/miFunctions.R", wdOA,wdOA_scripts))

muvTraits = c("FA_A","FA_B","SC_A","SC_B")

BioMetric_multivariate = BioMetric_Pheno%>%
  select(c(Sub,TwinType),muvTraits)

describe(BioMetric_multivariate[,muvTraits], skew=F)

# # Recode Data for Analysis - Rescale variables to have variances around 1.0
# BioMetric_multivariate = BioMetric_multivariate%>%
#   mutate(
#     mm2_z_FAres_1 = FA_A * 4,
#     mm2_z_FAres_2 = FA_B * 4,
#     mm2_z_SCres_1 = SC_A * 4,
#     mm2_z_SCres_2 =  SC_B * 4
#   )



describe(BioMetric_multivariate[,muvTraits], skew=F)

# Select Variables for Analysis
vars      <- c('FA','SC') 
nv        <- 2       # number of variables
ntv       <- nv*2    # number of total variables
selVars   <- paste(vars,c(rep("A",nv),rep("B",nv)),sep="_")


# Select Data for Analysis
mzData    <- subset(BioMetric_multivariate, TwinType==1, selVars)
dzData    <- subset(BioMetric_multivariate, TwinType==3, selVars)

# Generate Descriptive Statistics
round(colMeans(mzData,na.rm=TRUE),4)
round(colMeans(dzData,na.rm=TRUE),4)
round(cov(mzData,use="complete"),4)
round(cov(dzData,use="complete"),4)
round(cor(mzData,use="complete"),4)
round(cor(dzData,use="complete"),4)

# Set Starting Values 
svMe      <- c(0,0,0)                  # start value for means
svVa      <- .3                        # start value for Variance a
svVe      <- .6                       # start value for Variance  e
lbPa      <- .001                       # lower bound for path coefficient
# lbPaD     <- valDiagLU(lbPa,-10,NA,nv) # lower bound for diagonal, lower & upper elements of covariance matrix

# ------------------------------------------------------------------------------
# PREPARE MODEL

# ACE Model
# Create Algebra for expected Mean Matrices
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svVa,nv), label=labLower("VA",nv), name="VA" ) 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv), label=labLower("VC",nv), name="VC" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svVe,nv), label=labLower("VE",nv), name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+VC, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list(meanG, covA, covC, covE, covP )
modelMZ   <- mxModel( name="MZ", pars, covMZ, expCovMZ, dataMZ, expMZ, funML )
modelDZ   <- mxModel( name="DZ", pars, covDZ, expCovDZ, dataDZ, expDZ, funML )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Standardization
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

# Calculate genetic and environmental correlations
corA      <- mxAlgebra( expression=cov2cor(VA), name ="rA" )
corC      <- mxAlgebra( expression=cov2cor(VC), name ="rC" )
corE      <- mxAlgebra( expression=cov2cor(VE), name ="rE" )

# Create Algebra for Unstandardized and Standardized Variance Components
rowUS <- rep('US',nv)
colUS <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
estUS <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )

# Create Confidence Interval Objects
#create indicies for lower matricies and diagonal V and S components
row1 = seq(1,12,2)
row2 = seq(1,12,1)
ciACE <- mxCI( c("US[1,row1]","US[2,row2]") )
# Create Algebra for Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE'),each=nv)
estcorr     <- mxAlgebra( expression=cbind(rA,rC,rE), name="corr", dimnames=list(rowcorr,colcorr) )
## Create Confidence Interval Objects
ciCorr      <- mxCI(c("corr"))


# Build Model with Confidence Intervals
calc      <- list( matI, invSD, corA, corC, corE,estUS, ciACE,estcorr,ciCorr)
modelACE  <- mxModel( "mulACEvc", pars, modelMZ, modelDZ, multi, calc )

# ------------------------------------------------------------------------------
# RUN MODEL
# Run ACE Model
fitACE    <- mxRun( modelACE, intervals=T )
sumACE    <- summary( fitACE )
fitEsts(fitACE)

# See parameter specification (i.e varcov matricies)
parameterSpecifications(fitACE)

# Compare with Saturated Model
mxCompare( fitSAT, fitACE )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACE)
fitEsts(fitACE)

# ------------------------------------------------------------------------------
# RUN SUBMODELS

# Run AE model
modelAE   <- mxModel( fitACE, name="mulAEvc" )
modelAE   <- omxSetParameters( modelAE, labels=labLower("VC",nv), free=FALSE, values=0 )
fitAE     <- mxRun( modelAE, intervals=T )
mxCompare( fitACE, fitAE )
fitGofs(fitAE)

# Test if variances for shared E term can be removed
modelAE_1    <- mxModel( fitAE, name="mulAEvc_1" )
modelAE_1   <- omxSetParameters( modelAE_1, labels=c("VE21"), free=FALSE, values=0 )
fitAE_1      <- mxRun( modelAE_1, intervals=T )
mxCompare( fitAE, fitAE_1 )
# Test if variances for shared A term can be removed
modelAE_2    <- mxModel( fitAE, name="mulAEvc_2" )
modelAE_2   <- omxSetParameters( modelAE_1, labels=c("VA21"), free=FALSE, values=0 )
fitAE_2      <- mxRun( modelAE_2, intervals=T )
mxCompare( fitAE, fitAE_2 )


# Print Comparative Fit Statistics
mxCompare( fitACE, nested <- list(fitAE,fitAE_1,fitAE_2 ) )


#final model
fitGofs(fitAE)
round(rbind(fitACE$US$result,fitAE$US$result),4)
fitAE$US$result


# ACE Covariance Matrices & Proportions of Variance Matrices
matAEcov   <- c("VA","VE","V","VA/V","VE/V")
labAEcov   <- c("covA","covE","Var","stCovA","stCovE")
formatOutputMatrices(fitAE, matAEcov, labAEcov, vars,4)

# ACE Correlation Matrices 
matAEcor   <- c("solve(sqrt(I*VA)) %&% VA","solve(sqrt(I*VE)) %&% VE")
labAEcor   <- c("corA","corE")
formatOutputMatrices(fitAE, matAEcor, labAEcor, vars, 4)


#Store outputs####
#_h and bi_h estimates
AEvc = as.data.frame(round(fitAE$US$result,3))%>%
  as.data.frame()
colnames(AEvc) = c("VA_scenes","VA_faces",
                   "VC_scenes","VC_faces",
                   "VE_scenes","VE_faces",
                   "SA_scenes","SA_faces",
                   "SC_scenes","SC_faces",
                   "SE_scenes","SE_faces")
rownames(AEvc) = c("scenes","faces")
AEvc = AEvc%>%mutate(sample = "Sutherland et al. 2020",
                     facet = "e-b res.")
#_confidence interval
#create a matrix to store outputs
AECis = fitEstCis(fitAE)[1:18,]


#create newrownames (lower matrix)
namesCimatrix = c("VA11","VC11","VE11","SA11","SC11","SE11", #first row
                  "VA21","VA22","VC21","VC22","VE21","VE22","SA21","SA22","SC21","SC22","SE21","SE22" #second row
)

rownames(AECis) = namesCimatrix 
AECis = AECis%>%as.data.frame()%>%rownames_to_column()%>%mutate(sample = "Sutherland et al. 2020",
                                                                facet = "e-b res.")


#_genetic and environmental correlations
#confidence intervals
AECr_ci = fitEstCis(fitAE)[19:nrow(fitEstCis(fitAE)),]

#genetic correlations matrix
rA_mat = fitAE$rA$result%>%as.data.frame()%>%
  mutate(sample = "Sutherland et al. 2020",
         facet = "e-b res.",
         correlation = "rA")
colnames(rA_mat) = c("scenes","faces","sample","facet","correlation")
rownames(rA_mat) = c("scenes","faces")

#genetic correlations with confidence intervals
rA = rbind(
  AECr_ci%>%as.data.frame()%>%rowid_to_column()%>%filter(rowid == 2)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-faces")
)%>%
  mutate(sample = "Sutherland et al. 2020",
         facet = "e-b res.",
         correlation = "rA")

#environmental correlations matrix
rE_mat =fitAE$rE$result%>%as.data.frame()%>%
  mutate(sample = "Sutherland et al. 2020",
         facet = "e-b res.",
         correlation = "rE")
colnames(rE_mat) = c("scenes","faces","sample","facet","correlation")
rownames(rE_mat) = c("scenes","faces")

#environmental correlations with confidence intervals
rE = rbind(
  AECr_ci%>%as.data.frame()%>%rowid_to_column()%>%filter(rowid == 10)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-faces")
)%>%
  mutate(sample = "Sutherland et al. 2020",
         facet = "e-b res.",
         correlation = "rE")

mul_r = rbind(rA,rE)
####SAVE###
write_csv(AEvc,sprintf("%s/%s/02_Sutherland_2020/23_AE_ttres_multivariate_Sutherland2020_h2.csv",wdOA,wdOA_output))#save estimates
write_csv(AECis,sprintf("%s/%s/02_Sutherland_2020/23_AE_ttres_multivariate_Sutherland2020_CI.csv",wdOA,wdOA_output))#save ci
write_csv(mul_r,sprintf("%s/%s/02_Sutherland_2020/23_AE_ttres_multivariate_Sutherland2020_rArE.csv",wdOA,wdOA_output))#save correlations