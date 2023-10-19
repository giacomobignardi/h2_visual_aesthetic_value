#Author: Giacomo Bignardi
#Adapted from: Hermine Maes 01 04 2018
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description:
# Twin Multivariate (trivariate) CTD Structural Equation Model: 2 groups (MZss and DZss)
# Matrix style model - Raw data - Continuous data
#clean working enviroment 
rm(list = ls())
library(tidyverse)
library(OpenMx)
library(psych)
library(readr)
mxOption(NULL,"Default optimizer","SLSQP")

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"
OA_ImageOutput = "05_images/image"

#load dataFrames:
BioMetric_Pheno  = read_csv(sprintf("%s/%s/01_Germine_2015/05_twin_eb_multivariate.csv", wdOA,wdOA_output))
#load fit sat to compare with saturated model
load(sprintf("%s/%s/01_Germine_2015/12_1_SAT_eb_multivariate_Germine2015.RData", wdOA,wdOA_output))

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

# Select Data for Analysis
mzData    <- subset(BioMetric_multivariate, Zygosity==1, c(selVars)) #covVars
dzData    <- subset(BioMetric_multivariate, Zygosity==3, c(selVars)) #covVars

# Generate Descriptive Statistics
round(colMeans(mzData,na.rm=TRUE),4)
round(colMeans(dzData,na.rm=TRUE),4)
round(cov(mzData,use="complete"),4)
round(cov(dzData,use="complete"),4)
round(cor(mzData,use="complete"),4)
round(cor(dzData,use="complete"),4)

# Set Starting Values 
svBe      <- .0
svMe      <- c(0,0,0)                 # start value for means
svVa      <- .3                       # start value for Variance a
svVe      <- .4                       # start value for Variance  e
# ------------------------------------------------------------------------------
# PREPARE MODEL
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )
# expMean <- mxAlgebra( expression= meanG + cbind(defLage%*%blage,defLage%*%blage) + cbind(defLsex%*%blsex,defLsex%*%blsex), name="expMeanG" )
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
pars      <- list( meanG, covA, covC, covE, covP ) #pathBlage, pathBlsex,
modelMZ   <- mxModel( name="MZ", pars, covMZ, expCovMZ, dataMZ, expMZ, funML ) #defs
modelDZ   <- mxModel( name="DZ", pars, covDZ, expCovDZ, dataDZ, expDZ, funML ) #defs
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
row1 = seq(1,18,3)
row2 =  sort(c(seq(1,18,3),seq(2,18,3)))
row3 = seq(1,18,1)
ciACE <- mxCI( c("US[1,row1]","US[2,row2]","US[3,row3]") )

# Create Algebra for Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE'),each=nv)
estcorr     <- mxAlgebra( expression=cbind(rA,rC,rE), name="corr", dimnames=list(rowcorr,colcorr) )
## Create Confidence Interval Objects
ciCorr      <- mxCI(c("corr"))

# Build Model with Confidence Intervals
calc      <- list( matI, invSD, corA, corC, corE,estUS, ciACE ,estcorr,ciCorr)
modelACE  <- mxModel( "mulACEvc", pars, modelMZ, modelDZ, multi, calc )

# ------------------------------------------------------------------------------
# RUN MODEL
# Run ACE Model
fitACE    <- mxRun( modelACE, intervals=T )
sumACE    <- summary( fitACE )
fitEsts(fitACE)
fitEstCis(fitACE)

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
fitEsts(fitAE)

mxCompare( fitACE, fitAE )


# Run AE unique model
modelAE_Aunique   <- mxModel( fitAE, name="mulAEvcau" )
modelAE_Aunique   <- omxSetParameters( modelAE_Aunique, labels=c("VA21","VA31","VA32"), free=FALSE, values=0 )
fitAE_Aunique     <- mxRun( modelAE_Aunique, intervals=T )
mxCompare( fitAE, fitAE_Aunique )

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
#_model comparison
modelComparison = rbind(
  mxCompare( fitSAT, fitACE )%>% as.data.frame(),
  mxCompare( fitACE, nested <- list(fitAE, fitAE_Aunique) ) %>% as.data.frame() 
) %>% 
  mutate(sample = "Germine et al. 2015",
         facet = "e-b")%>% 
  filter(!row_number() %in% c(3))

#_h and bi_h estimates
AEvc = as.data.frame(round(fitAE$US$result,3))%>%
  as.data.frame()
colnames(AEvc) = c("VA_scenes","VA_faces","VA_abstracts",
                   "VC_scenes","VC_faces","VC_abstracts",
                   "VE_scenes","VE_faces","VE_abstracts",
                   "SA_scenes","SA_faces","SA_abstracts",
                   "SC_scenes","SC_faces","SC_abstracts",
                   "SE_scenes","SE_faces","SE_abstracts")
rownames(AEvc) = c("scenes","faces","abstracts")
AEvc = AEvc%>%mutate(sample = "Germine et al. 2015",
       facet = "e-b")
#_confidence interval
#create newrownames (lower matrix)
AECis = fitEstCis(fitAE) %>% as.data.frame()%>%rownames_to_column() %>% filter(grepl("US",rowname))
namesCimatrix = c("VA11",              "VC11",              "VE11",              "SA11",              "SC11",              "SE11",              #first row
                  "VA21","VA22",       "VC21","VC22",       "VE21","VE22",       "SA21","SA22",       "SC21","SC22",       "SE21","SE22",       #second row
                  "VA31","VA32","VA33","VC31","VC32","VC33","VE31","VE32","VE33","SA31","SA32","SA33","SC31","SC32","SC33","SE31","SE32","SE33" #third row
)
rownames(AECis) = namesCimatrix 
AECis = AECis%>%select(-rowname)%>%rownames_to_column()%>%mutate(sample = "Germine et al. 2015",
                      facet = "e-b")

#_genetic and environmental correlations
#confidence intervals
AECr_ci = fitEstCis(fitAE)%>% as.data.frame()%>%rownames_to_column() %>% filter(grepl("corr",rowname))

#genetic correlations matrix
rA_mat = fitAE$rA$result%>%as.data.frame()%>%
  mutate(sample = "Germine et al. 2015",
         facet = "e-b",
         correlation = "rA")
colnames(rA_mat) = c("scenes","faces","abstracts","sample","facet","correlation")
rownames(rA_mat) = c("scenes","faces","abstracts")

#genetic correlations with confidence intervals
rA = rbind(
  AECr_ci%>%rowid_to_column()%>%filter(rowid == 2)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-faces"),
  AECr_ci%>%rowid_to_column()%>%filter(rowid == 3)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-abstracts"),
  AECr_ci%>%rowid_to_column()%>%filter(rowid == 6)%>%rename(trait = "rowid")%>%mutate(trait ="abstracts-faces")
)%>%
  mutate(sample = "Germine et al. 2015",
         facet = "e-b",
         correlation = "rA")

#environmental correlations matrix
rE_mat =fitAE$rE$result%>%as.data.frame()%>%
  mutate(sample = "Germine et al. 2015",
         facet = "e-b",
         correlation = "rE")
colnames(rE_mat) = c("scenes","faces","abstracts","sample","facet","correlation")
rownames(rE_mat) = c("scenes","faces","abstracts")

#environmental correlations with confidence intervals
rE = rbind(
  AECr_ci%>%rowid_to_column()%>%filter(rowid == 20)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-faces"),
  AECr_ci%>%rowid_to_column()%>%filter(rowid == 21)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-abstracts"),
  AECr_ci%>%rowid_to_column()%>%filter(rowid == 24)%>%rename(trait = "rowid")%>%mutate(trait ="faces-abstracts")
)%>%
  mutate(sample = "Germine et al. 2015",
         facet = "e-b",
         correlation = "rE")

mul_r = rbind(rA,rE)

####SAVE###
save(fitAE, file = sprintf("%s/%s/01_Germine_2015/13_1.2_AE_eb_multivariate_Germine2015_SLSQP.RData",wdOA,wdOA_output))#save output
write_csv(modelComparison,sprintf("%s/%s/01_Germine_2015/13_1.2_AE_eb_multivariate_modelComparison_SLSQP.csv",wdOA,wdOA_output))#save model comparison
write_csv(AEvc,sprintf("%s/%s/01_Germine_2015/13_1.2_AE_eb_multivariate_h2_SLSQP.csv",wdOA,wdOA_output))#save estimates
write_csv(AECis,sprintf("%s/%s/01_Germine_2015/13_1.2_AE_eb_multivariate_CI_SLSQP.csv",wdOA,wdOA_output))#save ci
write_csv(mul_r,sprintf("%s/%s/01_Germine_2015/13_1.2_AE_eb_multivariate_rArE_SLSQP.csv",wdOA,wdOA_output))#save correlations
