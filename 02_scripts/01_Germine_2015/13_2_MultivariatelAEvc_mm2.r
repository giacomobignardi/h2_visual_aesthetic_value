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

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"
OA_ImageOutput = "05_images/image"

#load dataFrames:
BioMetric_Pheno  = read_csv(sprintf("%s/%s/01_Germine_2015/05_twin_tt_multivariate.csv", wdOA,wdOA_output))
#load fit sat to compare with saturated model
load(sprintf("%s/%s/01_Germine_2015/12_2_SAT_tt_multivariate_Germine2015.RData", wdOA,wdOA_output))

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
svMe      <- c(0,0,0)                 # start value for means
svVa      <- .5                       # start value for Variance a
svVe      <- .5                       # start value for Variance  e
# lbPa      <- .01                      # lower bound for path coefficient
# lbPaD     <- valDiagLU(lbPa,-10,NA,nv)# lower bound for diagonal, lower & upper elements of covariance matrix

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
#create indices for lower matrices and diagonal V and S components
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
calc      <- list( matI, invSD, corA, corC, corE, estUS, ciACE, estcorr,ciCorr )
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
# RUN EXPECTED SUB-MODEL
# Test if variances for A fo Abstract Objects term can be removed
modelCE_AO    <- mxModel( fitACE, name="mulCEao" )
modelCE_AO   <- omxSetParameters( modelCE_AO, labels=c("VA31","VA32","VA33"), free=FALSE, values=0 )
fitCE_AO     <- mxRun( modelCE_AO, intervals=T )
#Test if variances for C for scenes term can be removed
modelAE_SC    <- mxModel( fitCE_AO, name="mulAEsc" )
modelAE_SC   <- omxSetParameters( fitCE_AO, labels=c("VC11","VC21","VC31"), free=FALSE, values=0 )
fitAE_SC      <- mxRun( modelAE_SC, intervals=T )
#Test if variances for C for Faces term can be removed
modelAE_FA    <- mxModel( fitAE_SC, name="mulAEfa" )
modelAE_FA   <- omxSetParameters( modelAE_FA, labels=c("VC21","VC22","VC32"), free=FALSE, values=0 )
fitAE_FA      <- mxRun( modelAE_FA, intervals=T )
#Test if shared variances for A for natural images can be removed
modelAE_AOnoA    <- mxModel( fitAE_FA, name="mulAEscnoa" )
modelAE_AOnoA   <- omxSetParameters( modelAE_AOnoA, labels=c("VA21"), free=FALSE, values=0 )
fitAE_AOnoA    <- mxRun( modelAE_AOnoA, intervals=T )
#Test if shared variances for E for abstract objects can be removed
modelCE_AOnoE    <- mxModel( fitAE_FA, name="mulAEscnoe" )
modelCE_AOnoE   <- omxSetParameters( modelCE_AOnoE, labels=c("VE21","VE31"), free=FALSE, values=0 )
fitCE_AOnoE    <- mxRun( modelCE_AOnoE, intervals=T )

fitFinal = fitAE_FA
fitEsts(fitFinal)
fitEstCis(fitACE)

AECis = fitAE_FA$US$result
#ACE Covariance Matrices & Proportions of Variance Matrices
matAEcov   <- c("VA","VC","VE","V","VA/V","VC/V","VE/V")
labAEcov   <- c("covA","covC","covE","Var","stCovA","stCovC","stCovE")
formatOutputMatrices(fitFinal, matAEcov, labAEcov, vars,4)

#Store outputs####
#_model comparison
modelComparison =rbind(
  mxCompare( fitSAT, fitACE )%>% as.data.frame(),
  mxCompare( fitACE, nested <- list(fitCE_AO, fitAE_SC,fitAE_FA, fitAE_AOnoA,fitCE_AOnoE) ) %>% as.data.frame() 
) %>% 
  mutate(sample = "Germine et al. 2015",facet = "t-t")%>% 
  filter(!row_number() %in% c(3))

#_h and bi_h estimates
AEvc = as.data.frame(round(fitFinal$US$result,3))%>%
  as.data.frame()
colnames(AEvc) = c("VA_scenes","VA_faces","VA_abstracts",
                   "VC_scenes","VC_faces","VC_abstracts",
                   "VE_scenes","VE_faces","VE_abstracts",
                   "SA_scenes","SA_faces","SA_abstracts",
                   "SC_scenes","SC_faces","SC_abstracts",
                   "SE_scenes","SE_faces","SE_abstracts")
rownames(AEvc) = c("scenes","faces","abstracts")
AEvc = AEvc%>%mutate(sample = "Germine et al. 2015",
                     facet = "t-t")
#_confidence interval
#create newrownames (lower matrix)
AECis = fitEstCis(fitFinal)%>% as.data.frame()%>%rownames_to_column() %>% filter(grepl("US",rowname))
namesCimatrix = c("VA11",              "VC11",              "VE11",              "SA11",              "SC11",              "SE11",              #first row
                  "VA21","VA22",       "VC21","VC22",       "VE21","VE22",       "SA21","SA22",       "SC21","SC22",       "SE21","SE22",       #second row
                  "VA31","VA32","VA33","VC31","VC32","VC33","VE31","VE32","VE33","SA31","SA32","SA33","SC31","SC32","SC33","SE31","SE32","SE33" #third row
)
rownames(AECis) = namesCimatrix 
AECis = AECis%>%select(-rowname)%>%rownames_to_column()%>%mutate(sample = "Germine et al. 2015",facet = "t-t")                                                 
#_genetic and environmental correlations
#confidence intervals
AECr_ci = fitEstCis(fitFinal)%>% as.data.frame()%>%rownames_to_column() %>% filter(grepl("corr",rowname))

#genetic correlations matrix
rA_mat = fitFinal$rA$result%>%as.data.frame()%>%
  mutate(sample = "Germine et al. 2015",
         facet = "t-t",
         correlation = "rA")
colnames(rA_mat) = c("scenes","faces","abstracts","sample","facet","correlation")
rownames(rA_mat) = c("scenes","faces","abstracts")

#genetic correlations with confidence intervals
rA = rbind(
AECr_ci%>%rowid_to_column()%>%filter(rowid == 4)%>%rename(trait = "rowid")%>%mutate(trait ="scenes-faces")
)%>%
  mutate(sample = "Germine et al. 2015",
         facet = "t-t",
         correlation = "rA")
#environmental correlations matrix
rE_mat =fitFinal$rE$result%>%as.data.frame()%>%
  mutate(sample = "Germine et al. 2015",
         facet = "t-t",
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
         facet = "t-t",
         correlation = "rE")

mul_r = rbind(rA,rE)

####SAVE###
save(fitAE_FA, file = sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_Germine2015.RData",wdOA,wdOA_output))#save output
write_csv(modelComparison,sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_modelComparison.csv",wdOA,wdOA_output))#save model comparison
write_csv(AEvc,sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_h2.csv",wdOA,wdOA_output))#save estimates
write_csv(AECis,sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_CI.csv",wdOA,wdOA_output))#save ci
write_csv(mul_r,sprintf("%s/%s/01_Germine_2015/13_2_AE_tt_multivariate_rArE.csv",wdOA,wdOA_output))#save correlations

# ####Visualize###
# get_lower_tri<-function(cormat){
#   cormat[upper.tri(cormat)] <- NA
#   return(cormat)
# }
# loweTri = get_lower_tri(round(resutlsAE[,10:12],2))
# melted_AECis = reshape2::melt(loweTri)
# 
# melted_AECis[3,]$value = NA
# melted_AECis[6,]$value = NA
# melted_AECis[9,]$value = NA
# 
# melted_AECis$Var1 = c(rep("mm2-pl",3), rep("mm2-fa",3), rep("mm2-ao",3) )
# melted_AECis$Var2 =c (rep(c("mm2-pl","mm2-fa", "mm2-ao"), 3))
# 
# 
# p1 = ggplot(data = melted_AECis, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(color = "white") +
#   geom_text(aes(Var1,Var2, label = value), color = "white", 
#             size = 4)+
#   scale_fill_viridis_c(
#     limits = c(0, 1), 
#     na.value = 'white') +  
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank())
# 
# 
# get_lower_tri<-function(cormat){
#   cormat[upper.tri(cormat)] <- NA
#   return(cormat)
# }
# 
# loweTri = get_lower_tri(round(resutlsAE[,16:18],2))
# melted_AECis = reshape2::melt(loweTri)
# 
# melted_AECis$Var1 = c(rep("mm2-pl",3), rep("mm2-fa",3), rep("mm2-ao",3) )
# melted_AECis$Var2 =c (rep(c("mm2-pl","mm2-fa", "mm2-ao"), 3))
# 
# 
# p2 = ggplot(data = melted_AECis, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(color = "white") +
#   geom_text(aes(Var1,Var2, label = value), color = "white", 
#             size = 4)+
#   scale_fill_viridis_c(
#     limits = c(0, 1), 
#     na.value = 'white') +  
#   theme_minimal(base_size = 14)+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank())
# 
# pdf(sprintf("%s/%s/13_2_multivariate_mm2_genetic_Germine2015.pdf",wdNOA,wdNOA_ImageOutput),
#     width = 4,
#     height = 3)
# p1
# dev.off()
# 
# pdf(sprintf("%s/%s/13_2_multivariate_mm2_environmental_Germine2015.pdf",wdNOA,wdNOA_ImageOutput),
#     width = 4,
#     height = 3)
# p2
# dev.off()
# 
# # ####CHECK####
# # # Fit Common Pathway AE Model
# # # ------------------------------------------------------------------------------                                           
# # nl        <- 1
# # svFl      <- 1
# # svFa      <- .4
# # svFc      <- .2
# # svFe      <- .6
# # svfl      <- .6
# # 
# # 
# # # Recreate Algebra for expected Variance/Covariance Matrices in MZ & DZ twins (VC = 0)
# # covP      <- mxAlgebra( expression= VA+VE, name="V" )
# # covMZ     <- mxAlgebra( expression= VA, name="cMZ" )
# # covDZ     <- mxAlgebra( expression= 0.5%x%VA, name="cDZ" )
# # expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
# # expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
# # 
# # # Matrices Al, Cl, and El for variance components for latent phenotype
# # covAl     <- mxMatrix( type="Symm", nrow=nl, ncol=nl, free=TRUE, values=svFl, label=labLower("Al",nl), name="Al" ) 
# # # covCl     <- mxMatrix( type="Symm", nrow=nl, ncol=nl, free=TRUE, values=svFl, label=labLower("Cl",nl), name="Cl" )
# # covEl     <- mxMatrix( type="Symm", nrow=nl, ncol=nl, free=TRUE, values=svFl, label=labLower("El",nl), name="El" )
# # 
# # # Matrix and Algebra for constraint on variance of latent phenotype
# # covarLP   <- mxAlgebra( expression= Al+El, name="CovarLP" )
# # varLP     <- mxAlgebra( expression= diag2vec(CovarLP), name="VarLP" )
# # unit      <- mxMatrix( type="Unit", nrow=nl, ncol=1, name="Unit")
# # varLP1    <- mxConstraint( expression=VarLP == Unit, name="varLP1")
# 
# # Matrix f for factor loadings on latent phenotype
# pathFl    <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=TRUE, values=svfl, labels=labFull("fl",nv,nl), name="fl" )
# 
# # Matrices As, Cs, and Es for variance components for specific factors
# covAs     <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svFa, labels=labDiag("As",nv), name="As" )
# # covCs     <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svFc, labels=labDiag("Cs",nv), name="Cs" )
# covEs     <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svFe, labels=labDiag("Es",nv), name="Es" )
# 
# # Matrices A, C, and E compute variance components
# covA      <- mxAlgebra( expression=fl %&% Al + As, name="VA" )
# # covC      <- mxAlgebra( expression=fl %&% Cl + Cs, name="VC" )
# covE      <- mxAlgebra( expression=fl %&% El + Es, name="VE" )
# 
# # Create Model Objects for Multiple Groups
# pars      <- list(meanG, matI, invSD,
#                   covAl, covEl, covarLP, varLP, unit, pathFl,covAs, covEs, covA, covE, covP)
# modelMZ   <- mxModel( name="MZ", pars, covMZ, expCovMZ, dataMZ, expMZ, funML )
# modelDZ   <- mxModel( name="DZ", pars, covDZ, expCovDZ, dataDZ, expDZ, funML )
# multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
# 
# # Build & Run Model 
# #full common pathway:
# modelCP   <- mxModel( "mulCPc", pars, varLP1, modelMZ, modelDZ, multi )
# fitCP     <- mxRun(modelCP, intervals=F )
# fitGofs(fitCP)
# fitEsts(fitCP)
# 
# #compare AIC 
# summary(fitCP);summary(fitCP)
# mxCompare( fitAE, fitCP )
# 
# # Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices
# matCPpaths <- c("Al","El","iSD %*% fl","iSD %*% As","iSD %*% Es")
# labCPpaths <- c("stPathAl","stPathEl","stPathFl","stPathAs","stPathEs")
# formatOutputMatrices(fitCP, matCPpaths, labCPpaths, vars, 4)
# 
# 
# 
# # Fit Independent Pathway ACE Model
# # ------------------------------------------------------------------------------                                           
# nl        <- 2      # number of factors
# 
# # Matrices Al, El and fl for variance components for common factors (C set to 0)
# covAl     <- mxMatrix( type="Diag", nrow=nl, ncol=nl, free=c(F,F), values=c(1,0), label=labDiag("Al",nl), name="Al" ) 
# covEl     <- mxMatrix( type="Diag", nrow=nl, ncol=nl, free=c(F,F), values=c(0,1), label=labDiag("El",nl), name="El" )
# pathFl    <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=TRUE, values=svfl, labels= c("fla11","fla12","fla13","fle11","fle12","fle13"), name="fl" )
# 
# # Create Model Objects for Multiple Groups
# pars      <- list(meanG, matI, invSD,
#                   covAl, covEl, covarLP, pathFl,covAs, covEs, covA, covE, covP)
# modelMZ   <- mxModel( name="MZ", pars, covMZ, expCovMZ, dataMZ, expMZ, funML )
# modelDZ   <- mxModel( name="DZ", pars, covDZ, expCovDZ, dataDZ, expDZ, funML )
# multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
# 
# # Build & Run Model 
# modelIP   <- mxModel( "mulIPc", pars, modelMZ, modelDZ, multi )
# fitIP     <- mxRun(modelIP, intervals=F)
# 
# fitGofs(fitIP)
# fitEsts(fitIP)
# 
# #Compare
# summary(fitIP);summary(fitAE)
# 
# parameterSpecifications(fitIP)
# 
# # Generate List
# formatOutputMatrices(fitIP, matCPpaths, labCPpaths, vars, 4)

