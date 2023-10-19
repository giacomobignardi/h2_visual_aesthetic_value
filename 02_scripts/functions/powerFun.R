# Functions to do power analysis for Twins Models
#
#
#
# acePow(add, com, Nmz, Ndz)                                          # Univariate Continuous power analysis
# bivPow(A1, A2, rg, C1, C2, rc, Nmz, Ndz)                            # Bivariate Continuous power analysis
# powerPlot(maxN, Wncp, wu = T)                                       # Function to plot the power curves and sample sizes for 20% , 50% & 80% power. 
# powerValue(N, Wncp, wu = T)                                         # Function to calculate the power for a specific sample size
# props(prev,r)                                                       # Function to calculate the proportions in each cell for binary data
# OrdData(r, percent, N)                                              # Function to simulate univariate ordinal data
# BivOrdData(r, percent1, percent2, N)                                # Function to simulate bivariate ordinal data
# acePowOrdfunction(add, com, percent, Nmz, Ndz)                      # Function to calculate power for univariate ordinal twin models
# aceBivPowOrd(A1, A2, rg, C1, C2, rc, percent1, percent2, Nmz, Ndz)  # Function to calculate power for bivariate ordinal twin models
# sexLimPower(am, cm, af, cf, rg, Nmzm, Nmzf, Ndzm, Ndzf, Ndzo)       # Funciton to calculate power for continuous sex limitation models
# CholSim(aLower, cLower, eLower, Nmz, Ndz)                           # Function to simulate continuous Cholesky decomposition data
# indPathSim(aLoad, cLoad, eLoad, specA, specC, specE, Nmz, Ndz)      # Function to simulate Independent Pathway Twin Data
# comPathSim(aL, cL, eL, load, specA, specC, specE, Nmz, Ndz)         # Function to simulate Common Pathway Twin Data




acePow <- function(add, com, Nmz, Ndz){

nv <- 1
ntv <- nv*2

AA <- add
CC <- com
EE <- 1 - add - com

mzMat <- matrix(c(AA + CC + EE, AA + CC, AA + CC, AA + CC + EE),2)
dzMat <- matrix(c(AA + CC + EE, .5*AA + CC, .5*AA + CC, AA + CC + EE),2)

mzData <- mvrnorm(Nmz , mu = c(0,0), mzMat, empirical = T)
dzData <- mvrnorm(Ndz , mu = c(0,0), dzMat, empirical = T)

selVars <- paste("t", 1:2, sep = "")
colnames(mzData) <- colnames(dzData) <- selVars

MZdata  <-  mxData(mzData, type="raw" )
DZdata  <-  mxData(dzData, type="raw" )


A <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.6, label="A11", name="A" )
C <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.6, label="C11", name="C" )
E <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.6, label="E11", name="E" )

Mean    <-  mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= 0, label="mean", name="Mean" )
expMean <-  mxAlgebra( expression= cbind(Mean,Mean), name="expMean")

expCovMZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),cbind(A+C   , A+C+E)), name="expCovMZ" )
expCovDZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 

obs <-  list(A,C,E,Mean,expMean,expCovMZ,expCovDZ)

fun <- mxFitFunctionML()
mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars )
dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars )

MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)

aceFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
ace <- mxModel("ACE", MZ, DZ, fun, aceFun)

aceFit <- suppressWarnings(mxRun(ace, silent = T))
summary(aceFit)

ceFit <- omxSetParameters(aceFit, labels = "A11", free = F, values = 0)
ceFit <- suppressWarnings(mxRun(ceFit, silent = T))

aeFit <- omxSetParameters(aceFit, labels = "C11", free = F, values = 0)
aeFit <- suppressWarnings(mxRun(aeFit, silent = T))

ncpA <- (ceFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
ncpC <- (aeFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])


Ests <-    as.data.frame(cbind(rbind(add, com), rbind(aceFit$output$estimate[1], aceFit$output$estimate[2])))
colnames(Ests) <- c("Simulated", "Estimated")
Ests$Difference <- Ests$Simulated - Ests$Estimated
Ests <- round(Ests, 3)
Imprecision <- sum(Ests$Difference)


if( Imprecision > .02 ) warning("The estimated values may differ from the simulated values more than you are comfortable with. 
Make sure you double check these values. You may want to increase the sample size ", call. = F)
return ( list(Parameters = Ests, WncpA = ncpA, WncpC =ncpC))

}



bivPow <- function(A1, A2, rg, C1, C2, rc, Nmz, Ndz){

nv <- 2
ntv <- nv*2

aCOR <- matrix(c(1, rg, rg, 1),2)
Avar <- c(A1, A2)
cCOR <- matrix(c(1, rc, rc, 1),2)
Cvar <- c(C1, C2)


AA <- diag(sqrt(Avar)) %&% aCOR
CC <- diag(sqrt(Cvar)) %&% cCOR
EE <- vec2diag(diag(diag(2)) - diag(AA) - diag(CC))

mzMat <- rbind(cbind(AA + CC + EE, AA + CC), cbind (AA + CC, AA + CC + EE))
dzMat <- rbind(cbind(AA + CC + EE, .5* AA + CC), cbind (.5*AA + CC, AA + CC + EE))

mzData <- mvrnorm(Nmz , mu = c(0,0,0,0), mzMat, empirical = T)
dzData <- mvrnorm(Ndz , mu = c(0,0,0,0), dzMat, empirical = T)

selVars <- paste(rep(c("t1", "v2"), each = 2), 1:2, sep = "_")
colnames(mzData) <- colnames(dzData) <- selVars

MZdata  <-  mxData(mzData, type="raw" )
DZdata  <-  mxData(dzData, type="raw" )



A <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=c(.33, 0, .33), label=c("A11","A21","A22"), name="A" )
C <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=c(.33, 0, .33), label=c("C11","C21","C22"), name="C" )
E <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=c(.34, 0, .34), label=c("E11","E21","E22"), name="E" )
expMean    <-  mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values= 0, label=c("mean1","mean2"), name="expMean" )
expCovMZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C), cbind(A+C   , A+C+E)), name="expCovMZ" )
expCovDZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C), cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 
obs <-  list(A,C,E,expMean,expCovMZ,expCovDZ)
fun <- mxFitFunctionML()
mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars )
dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars )

MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)

aceFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
ace <- mxModel("ACE", MZ, DZ, fun, aceFun)

aceFit <- suppressWarnings(mxRun(ace, silent = T))


ceFit1 <- omxSetParameters(aceFit, labels = "A11", free = F, values = 0)
ceFit1 <- suppressWarnings(mxRun(ceFit1, silent = T))

ceFit2 <- omxSetParameters(aceFit, labels = "A22", free = F, values = 0)
ceFit2 <- suppressWarnings(mxRun(ceFit2, silent = T))

rgFit <- omxSetParameters(aceFit, labels = "A21", free = F, values = 0)
rgFit <- suppressWarnings(mxRun(rgFit, silent = T))

aeFit1 <- omxSetParameters(aceFit, labels = "C11", free = F, values = 0)
aeFit1 <- suppressWarnings(mxRun(aeFit1, silent = T))

aeFit2 <- omxSetParameters(aceFit, labels = "C22", free = F, values = 0)
aeFit2 <- suppressWarnings(mxRun(aeFit2, silent = T))

rcFit <- omxSetParameters(aceFit, labels = "C21", free = F, values = 0)
rcFit <- suppressWarnings(mxRun(rcFit, silent = T))

ncpA.1  <- (ceFit1$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
ncpA.2  <- (ceFit2$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
ncpRg.1 <- (rgFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])

ncpC.1  <- (aeFit1$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
ncpC.2  <- (aeFit2$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
ncpRc.1 <- (rcFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])


pars <- round(cbind(aceFit$output$estimate, aceFit$output$standardErrors), 3)
colnames(pars) <- c("Estimates", "Std. Err.")

matA <- matrix(c(aceFit$output$estimate[1], aceFit$output$estimate[2],aceFit$output$estimate[2],aceFit$output$estimate[3]),2)
matC <- matrix(c(aceFit$output$estimate[4], aceFit$output$estimate[5],aceFit$output$estimate[5],aceFit$output$estimate[6]),2)

RG <- (solve(sqrt(diag(2) * matA)) %&% matA)[2,1]
RC <- (solve(sqrt(diag(2) * matC)) %&% matC)[2,1]

Ests <-    as.data.frame(cbind(	rbind(A1, A2, rg, C1, C2, rc), 
                                rbind(aceFit$output$estimate[1],aceFit$output$estimate[3], RG, aceFit$output$estimate[4], aceFit$output$estimate[6], RC)))
colnames(Ests) <- c("Simulated", "Estimated")
Ests$Difference <- Ests$Simulated - Ests$Estimated
Ests <- round(Ests, 3)
Imprecision <- sum(Ests$Difference)


if( Imprecision > .02 ) warning("The estimated values may differ from the simulated values more than you are comfortable with. 
Make sure you double check these values. You may want to increase the sample size ", call. = F)
return ( list(Parameters = Ests, WncpA1 = ncpA.1, WncpA2 = ncpA.2, WncpRg = ncpRg.1,
	WncpC1 = ncpC.1, WncpC2 = ncpC.2, WncpRc = ncpRc.1))
}



powerPlot <- function(maxN, Wncp, wu = T){
	n<- 1:maxN
	ncp <- matrix(NA, maxN, length(Wncp))
	pcrit <- ifelse(wu, .10, .05)
	for(i in 1:length(Wncp)){
	ncp[,i] <- 	Wncp[i] * n
	}	
			
plot(1- pchisq(qchisq(1- pcrit, 1), 1, ncp[,1]), type = "l", lwd = 3, col = 1, main = " ", ylab = "Power to Detect Significant Parameter" , xlab = "Sample Size",  ylim = c(0, 1), xlim = c(0, maxN))

for(j in 2:length(Wncp)){
	lines(1- pchisq(qchisq(1- pcrit, 1), 1, ncp[,j]), lwd = 3, col = j)
	}
	
	PowerTable <- matrix(NA, 3, dim(ncp)[2])
	for(k in 1:dim(ncp)[2]){
	PowerTable[,k] <- c(sum(1- pchisq(qchisq(1- pcrit, 1), 1, ncp[,k])<.2) + 1, 
	                    sum(1- pchisq(qchisq(1- pcrit, 1), 1, ncp[,k])<.5) + 1, 
						sum(1- pchisq(qchisq(1- pcrit, 1), 1, ncp[,k])<.8) + 1)
	}
	rownames(PowerTable ) <- c(".2", ".5", ".8")	
	colnames(PowerTable ) <- paste("N", 1:dim(ncp)[2], sep = "")	
	
	PowerTable[PowerTable == (maxN+1)] <- paste("More than", maxN)
	
	abline(v = PowerTable[3,], h = .8, lwd = .5, lty = 2, col = "gray")
	
	PowerTable
}


powerValue <- function(N, Wncp, wu = T){
	pcrit <- ifelse(wu, .10, .05)
	ncp <- Wncp * N
	pow <- 1- pchisq(qchisq(1- pcrit, 1), 1, ncp)
	pow
}

#######




props <- function(prev,r){
	thr <- qnorm(1-prev)
prop <-	matrix(	omxAllInt(
					matrix(c(1,r,r,1),2,2), 
					matrix(0,1,2), 
					t(matrix(c(rep(-Inf,2),thr,thr,rep(Inf,2)),2,3))), 2,2)					
prop				
}

OrdData <- function(r, percent, N){
	if(sum(percent)!=1) stop("The cell percentages for Variable 1 do not sum to 1.")

thr <- c(-Inf, qnorm(cumsum(percent)))

props<-round(omxAllInt(matrix(c(1,r,r,1),2,2), 
		  matrix(0,1,2), 
		  matrix(c(thr,thr),length(thr)))*N
		  ) 

combs	<- expand.grid(t1_x1 = 1:length(percent), t2_x1 = 1:length(percent))
combs <- combs[order(combs[,2], combs[,1]),]	

for(i in 1:length(props)){
obs<-	(matrix(combs[i,], round(props[i,]),2, byrow = T))
ifelse(i==1, dat <- obs, dat<- rbind(dat, obs))
}	

colnames(dat) <- colnames(combs)

dat				
	}

BivOrdData <- function(r, percent1, percent2, N){
	if(sum(percent1)!=1) stop("The cell percentages for Variable 1 do not sum to 1.")
	if(sum(percent2)!=1) stop("The cell percentages for Variable 2 do not sum to 1.")

thr1 <- c(-Inf, qnorm(cumsum(percent1)))
thr2 <- c(-Inf, qnorm(cumsum(percent2)))
while(length(thr1)<length(thr2)) thr1 <- c(thr1, NA)
while(length(thr1)>length(thr2)) thr2 <- c(thr2, NA)


props<-round(omxAllInt(matrix(r,4,4), 
		  matrix(0,1,4), 
		  matrix(c(thr1,thr2,thr1,thr2),max(length(thr1),length(thr2)),4)) *N)

combs	<- expand.grid(t1_x1 = 1:length(percent1), t1_x2 = 1:length(percent2),t2_x1 = 1:length(percent1), t2_x2 = 1:length(percent2))
combs <- combs[order(combs[,1], combs[,2], combs[,3]),]	

for(i in 1:length(props)){
obs<-	suppressWarnings(matrix(combs[i,], round(props[i,]),4, byrow = T))
ifelse(i==1, dat <- obs, dat<- rbind(dat, obs))
}	

colnames(dat) <- colnames(combs)

dat				
	}



	acePowOrd <- function(add, com, percent, Nmz, Ndz){

	nv <- 1
	ntv <- nv*2

	AA <- add
	CC <- com
	EE <- 1 - add - com

	mzMat <- matrix(c(AA + CC + EE, AA + CC, AA + CC, AA + CC + EE),2)
	dzMat <- matrix(c(AA + CC + EE, .5*AA + CC, .5*AA + CC, AA + CC + EE),2)

mzData <- OrdData(mzMat[2,1], percent = percent, N = Nmz)
dzData <- OrdData(dzMat[2,1], percent = percent, N = Ndz)

nTH <- max(unlist(mzData), unlist(dzData))-1

mzData <- mxFactor(as.data.frame(mzData), levels = 1:(nTH+1))
dzData <- mxFactor(as.data.frame(dzData), levels = 1:(nTH+1))

	selVars <- paste("t", 1:2, sep = "")
	colnames(mzData) <- colnames(dzData) <- selVars

	MZdata  <-  mxData(mzData, type="raw" )
	DZdata  <-  mxData(dzData, type="raw" )


	A <-  mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, label="A11", name="A" )
	C <-  mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, label="C11", name="C" )
	E <-  mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.4, label="E11", name="E" )

	vars = mxAlgebra( A+C+E, name="vars" )
	con1 <- mxConstraint(vars ==1, name = "con1")

	Mean    <-  mxMatrix( type="Full", nrow=1, ncol=nv, free=F, values= 0, label="mean", name="Mean" )
	expMean <-  mxAlgebra( expression= cbind(Mean,Mean), name="expMean")

	th1 <- (1:nTH)-1
	ThrVals <- matrix(c(th1), nTH, 2)

	Ups   <- mxMatrix("Full", nTH, 2, free = !is.na(ThrVals), values = ThrVals, labels =  paste("th", 1:nTH, sep = ""), name = "Ups")

	expCovMZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),cbind(A+C   , A+C+E)), name="expCovMZ" )
	expCovDZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 

	obs <-  list(A,C,E,vars,con1,Mean,expMean,Ups,expCovMZ,expCovDZ)

	fun <- mxFitFunctionML()
	mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", thresholds = "Ups", dimnames=selVars )
	dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", thresholds = "Ups", dimnames=selVars )

	MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
	DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)

	aceFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
	ace <- mxModel("ACE", MZ, DZ, fun, aceFun)

	aceFit <- suppressWarnings(mxRun(ace, silent = T))
	summary(aceFit)
	
	ceFit <- omxSetParameters(aceFit, labels = "A11", free = F, values = 0)
	ceFit <- suppressWarnings(mxRun(ceFit, silent = T))

	aeFit <- omxSetParameters(aceFit, labels = "C11", free = F, values = 0)
	aeFit <- suppressWarnings(mxRun(aeFit, silent = T))

	ncpA  <- (ceFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
	ncpC  <- (aeFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])


	Ests <-    as.data.frame(cbind(rbind(add, com), rbind(aceFit$output$estimate[1], aceFit$output$estimate[2])))
	colnames(Ests) <- c("Simulated", "Estimated")
	Ests$Difference <- Ests$Simulated - Ests$Estimated
	Ests <- round(Ests, 3)
	Imprecision <- sum(Ests$Difference)

	if( Imprecision > .02 ) warning("The estimated values may differ from the simulated values more than you are comfortable with. 
	Make sure you double check these values. You may want to increase the sample size ", call. = F)
	return ( list(Parameters = Ests, WncpA = ncpA, WncpC = ncpC))
	
		}
	



	aceBivPowOrd <- function(A1, A2, rg, C1, C2, rc, percent1, percent2, Nmz, Ndz){

	nv <- 2
	ntv <- nv*2

	aCOR <- matrix(c(1, rg, rg, 1),2)
	Avar <- c(A1, A2)
	cCOR <- matrix(c(1, rc, rc, 1),2)
	Cvar <- c(C1, C2)


	AA <- diag(sqrt(Avar)) %&% aCOR
	CC <- diag(sqrt(Cvar)) %&% cCOR
	EE <- vec2diag(diag(diag(2)) - diag(AA) - diag(CC))

	mzMat <- rbind(cbind(AA + CC + EE, AA + CC), cbind (AA + CC, AA + CC + EE))
	dzMat <- rbind(cbind(AA + CC + EE, .5* AA + CC), cbind (.5*AA + CC, AA + CC + EE))


	datMZ <- BivOrdData(r=mzMat, percent1 = percent1, percent2 = percent2, N = Nmz)
	datDZ <- BivOrdData(r=dzMat, percent1 = percent1, percent2 = percent2, N = Ndz)

	nTH <- max(unlist(datMZ), unlist(datDZ))-1
			
	mzData <- mxFactor(as.data.frame(datMZ), levels = 1:(nTH+1))
	dzData <- mxFactor(as.data.frame(datDZ), levels = 1:(nTH+1))

	for (ord in 1:4) {	
	mzData[,ord] <- mxFactor((mzData[,ord]),levels=c(min(unlist(mzData[,ord]), na.rm = T )):c(max(unlist(mzData[,ord]), na.rm = T )))			
	dzData[,ord] <- mxFactor((dzData[,ord]),levels=c(min(unlist(dzData[,ord]), na.rm = T )):c(max(unlist(dzData[,ord]), na.rm = T )))			
				}


	selVars <- colnames(datMZ)
	colnames(mzData) <- colnames(dzData) <- selVars


	MZdata  <-  mxData(mzData, type="raw" )
	DZdata  <-  mxData(dzData, type="raw" )


	####################
	####################

	A <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=c(.33, 0, .33), label=c("A11","A21","A22"), name="A" )
	C <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=c(.33, 0, .33), label=c("C11","C21","C22"), name="C" )
	E <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=c(.34, 0, .34), label=c("E11","E21","E22"), name="E" )

	expMean    <-  mxMatrix( type="Full", nrow=1, ncol=ntv, free=F, values= 0, label=c("mean1","mean2"), name="expMean" )

	th1 <- 1:(max(unlist(datMZ[,1]))-1) -1
	th2 <- 1:(max(unlist(datMZ[,2]))-1) -1
	while(length(th1)<length(th2)) th1 <- c(th1, NA)
	while(length(th1)>length(th2)) th2 <- c(th2, NA)
	ThrVals <- matrix(c(th1, th2), nTH, 4)


	Ups   <- mxMatrix("Full", nTH, 4, free = !is.na(ThrVals), values = ThrVals, labels = paste(rep(paste("v", 1:nv, sep=""),each = nTH), paste("th", 1:nTH, sep = ""),sep = "_"), name = "Ups")


	expCovMZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C), cbind(A+C   , A+C+E)), name="expCovMZ" )
	expCovDZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C), cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 

	vars = mxAlgebra( A+C+E, name="vars" )
	con1 <- mxConstraint(vars[1,1] ==1, name = "con1")
	con2 <- mxConstraint(vars[2,2] ==1, name = "con2")


	obs <-  list(A,C,E,expMean, Ups, expCovMZ,expCovDZ, vars, con1, con2)
	fun <- mxFitFunctionML()
	mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", thresholds = "Ups", dimnames=selVars )
	dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", thresholds = "Ups", dimnames=selVars )



	MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
	DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)

	aceFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
	ace <- mxModel("ACE", MZ, DZ, fun, aceFun)

	aceFit <- suppressWarnings(mxRun(ace, silent = T))
	#summary(aceFit)


	ceFit1 <- omxSetParameters(aceFit, labels = "A11", free = F, values = 0)
	ceFit1 <- suppressWarnings(mxRun(ceFit1, silent = T))

	ceFit2 <- omxSetParameters(aceFit, labels = "A22", free = F, values = 0)
	ceFit2 <- suppressWarnings(mxRun(ceFit2, silent = T))

	rgFit <- omxSetParameters(aceFit, labels = "A21", free = F, values = 0)
	rgFit <- suppressWarnings(mxRun(rgFit, silent = T))

	aeFit1 <- omxSetParameters(aceFit, labels = "C11", free = F, values = 0)
	aeFit1 <- suppressWarnings(mxRun(aeFit1, silent = T))

	aeFit2 <- omxSetParameters(aceFit, labels = "C22", free = F, values = 0)
	aeFit2 <- suppressWarnings(mxRun(aeFit2, silent = T))

	rcFit <- omxSetParameters(aceFit, labels = "C21", free = F, values = 0)
	rcFit <- suppressWarnings(mxRun(rcFit, silent = T))

	ncpA.1  <- (ceFit1$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
	ncpA.2  <- (ceFit2$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
	ncpRg.1 <- (rgFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])

	ncpC.1  <- (aeFit1$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
	ncpC.2  <- (aeFit2$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])
	ncpRc.1 <- (rcFit$output$Minus2LogLikelihood - aceFit$output$Minus2LogLikelihood)/(dim(mzData)[1] + dim(dzData)[1])

	matA <- matrix(c(aceFit$output$estimate[1], aceFit$output$estimate[2],aceFit$output$estimate[2],aceFit$output$estimate[3]),2)
	matC <- matrix(c(aceFit$output$estimate[4], aceFit$output$estimate[5],aceFit$output$estimate[5],aceFit$output$estimate[6]),2)

	RG <- (solve(sqrt(diag(2) * matA)) %&% matA)[2,1]
	RC <- (solve(sqrt(diag(2) * matC)) %&% matC)[2,1]

	Ests <-    as.data.frame(cbind(	rbind(A1, A2, rg, C1, C2, rc), 
	                                rbind(aceFit$output$estimate[1],aceFit$output$estimate[3], RG, aceFit$output$estimate[4], aceFit$output$estimate[6], RC)))
	colnames(Ests) <- c("Simulated", "Estimated")
	Ests$Difference <- Ests$Simulated - Ests$Estimated
	Ests <- round(Ests, 3)
	Imprecision <- sum(Ests$Difference)

	pars <- round(cbind(aceFit$output$estimate), 3)
	if( Imprecision > .04 ) warning("The estimated values may differ from the simulated values more than you are comfortable with. 
	Make sure you double check these values. You may want to increase the sample size ", call. = F)
	return ( list(Parameters = pars, WncpA1 = ncpA.1, WncpA2 = ncpA.2, WncpRg = ncpRg.1,
		WncpC1 = ncpC.1, WncpC2 = ncpC.2, WncpRc = ncpRc.1))
	
		}
	
		
## Sex Limitation Model Power Function

sexLimPower <- function(am, cm, af, cf, rg, Nmzm, Nmzf, Ndzm, Ndzf, Ndzo)	{
	TotalN <- sum(Nmzm, Nmzf, Ndzm, Ndzf, Ndzo)
	AM <- sqrt(am)
	CM <- sqrt(cm)
	EM <- sqrt(1- AM^2- CM^2)

	AF <- sqrt(af)
	CF <- sqrt(cf)
	EF <- sqrt(1- AF^2 - CF^2)

	RG <-  rg

	mzmCov <- matrix(c(AM^2+CM^2+EM^2, AM^2+CM^2, AM^2+CM^2, AM^2+CM^2+EM^2),2,2)
	mzfCov <- matrix(c(AF^2+CF^2+EF^2, AF^2+CF^2, AF^2+CF^2, AF^2+CF^2+EF^2),2,2)

	dzmCov <- matrix(c(AM^2+CM^2+EM^2, .5*AM^2+CM^2, .5*AM^2+CM^2, AM^2+CM^2+EM^2),2,2)
	dzfCov <- matrix(c(AF^2+CF^2+EF^2, .5*AF^2+CF^2, .5*AF^2+CF^2, AF^2+CF^2+EF^2),2,2)

	dzoCov <- matrix(c(AM^2+CM^2+EM^2 , .5*RG*AM*AF + CM*CF, .5*RG*AM*AF + CM*CF, AF^2+CF^2+EF^2),2,2)


	mzmData <- mvrnorm(Nmzm, c(0,0), mzmCov, empirical = T)
	mzfData <- mvrnorm(Nmzf, c(0,0), mzfCov, empirical = T)
	dzmData <- mvrnorm(Ndzm, c(0,0), dzmCov, empirical = T)
	dzfData <- mvrnorm(Ndzf, c(0,0), dzfCov, empirical = T)
	dzoData <- mvrnorm(Ndzo, c(0,0), dzoCov, empirical = T)

	selVars <- colnames(mzmData) <- colnames(mzfData) <- colnames(dzmData) <- colnames(dzfData) <- colnames(dzoData) <- paste("T", 1:2, sep = "")



	# ------------------------------------------------------------------------------
	nv <- 1
	ntv <- 2
	# General Sex Limitation ACE Model

	# Matrices declared to store a, c, and e Path Coefficients
	pathAf    <- mxMatrix( "Lower", nrow=nv, ncol=nv, free=TRUE, values=AF, label="af11", lbound = .0001, name="af" ) 
	pathCf    <- mxMatrix( "Lower", nrow=nv, ncol=nv, free=TRUE, values=CF, label="cf11", lbound = .0001, name="cf" )
	pathEf    <- mxMatrix( "Lower", nrow=nv, ncol=nv, free=TRUE, values=EF, label="ef11", lbound = .0001, name="ef" )
	pathAm    <- mxMatrix( "Lower", nrow=nv, ncol=nv, free=TRUE, values=AM, label="am11", lbound = .0001, name="am" ) 
	pathCm    <- mxMatrix( "Lower", nrow=nv, ncol=nv, free=TRUE, values=CM, label="cm11", lbound = .0001, name="cm" )
	pathEm    <- mxMatrix( "Lower", nrow=nv, ncol=nv, free=TRUE, values=EM, label="em11", lbound = .0001, name="em" )
	pathRg    <- mxMatrix( "Lower", nrow=1, ncol=1, free=TRUE, values=1,    label="rg11", name="rg", ubound=.9999, lbound=0 )
		
	# Matrices generated to hold A, C, and E computed Variance Components
	covAf     <- mxAlgebra( af %*% t(af), name="Af" )
	covCf     <- mxAlgebra( cf %*% t(cf), name="Cf" ) 
	covEf     <- mxAlgebra( ef %*% t(ef), name="Ef" )
	covAm     <- mxAlgebra( am %*% t(am), name="Am" )
	covCm     <- mxAlgebra( cm %*% t(cm), name="Cm" ) 
	covEm     <- mxAlgebra( em %*% t(em), name="Em" )

	# Algebra to compute total variances and standard deviations (diagonal only)
	covPf     <- mxAlgebra( Af+Cf+Ef, name="Vf" )
	covPm     <- mxAlgebra( Am+Cm+Em, name="Vm" )

	# Algebras generated to hold Parameter Estimates and Derived Variance Components
	colVarsZf <- c('Af','Cf','Ef','SAf','SCf','SEf')
	estVarsZf <- mxAlgebra( cbind(Af,Cf,Ef,Af/Vf,Cf/Vf,Ef/Vf), name="VarsZf", dimnames=list(NULL,colVarsZf))
	colVarsZm <- c('Am','Cm','Em','SAm','SCm','SEm')
	estVarsZm <- mxAlgebra( cbind(Am,Cm,Em,Am/Vm,Cm/Vm,Em/Vm), name="VarsZm", dimnames=list(NULL,colVarsZm))

	# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
	meanGf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=0, label="meanf", name="expMeanGf" )
	meanGm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=0, label="meanm", name="expMeanGm" )
	meanGfm   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=0, label=c("meanf","meanm"), name="expMeanGfm" )
	covMZf    <- mxAlgebra( expression= rbind( cbind(Vf, Af+Cf), cbind(Af+Cf, Vf)), name="expCovMZf" )
	covDZf    <- mxAlgebra( expression= rbind( cbind(Vf, 0.5%x%Af+Cf), cbind(0.5%x%Af+Cf, Vf)), name="expCovDZf" )
	covMZm    <- mxAlgebra( expression= rbind( cbind(Vm, Am+Cm), cbind(Am+Cm, Vm)), name="expCovMZm" )
	covDZm    <- mxAlgebra( expression= rbind( cbind(Vm, 0.5%x%Am+Cm), cbind(0.5%x%Am+Cm, Vm)), name="expCovDZm" )
	CVfm      <- mxAlgebra( expression= 0.5%*%rg%x%(af%*%t(am))+cf%*%t(cm), name="CVfm" )
	CVmf      <- mxAlgebra( expression= 0.5%*%rg%x%(am%*%t(af))+cm%*%t(cf), name="CVmf" )
	covDZo    <- mxAlgebra( expression= rbind( cbind(Vf, CVfm), cbind(CVmf, Vm)), name="expCovDZo" )

	# Data objects for Multiple Groups
	dataMZf   <- mxData( observed=mzfData, type="raw" )
	dataDZf   <- mxData( observed=dzfData, type="raw" )
	dataMZm   <- mxData( observed=mzmData, type="raw" )
	dataDZm   <- mxData( observed=dzmData, type="raw" )
	dataDZo   <- mxData( observed=dzoData, type="raw" )

	# Expectation objects for Multiple Groups
	expMZf    <- mxExpectationNormal( covariance="expCovMZf", means="expMeanGf", dimnames=selVars )
	expDZf    <- mxExpectationNormal( covariance="expCovDZf", means="expMeanGf", dimnames=selVars )
	expMZm    <- mxExpectationNormal( covariance="expCovMZm", means="expMeanGm", dimnames=selVars )
	expDZm    <- mxExpectationNormal( covariance="expCovDZm", means="expMeanGm", dimnames=selVars )
	expDZo    <- mxExpectationNormal( covariance="expCovDZo", means="expMeanGfm", dimnames=selVars )
	funML     <- mxFitFunctionML()

	# Combine Groups
	parsZf    <- list( pathAf, pathCf, pathEf, covAf, covCf, covEf, covPf, estVarsZf )
	parsZm    <- list( pathAm, pathCm, pathEm, covAm, covCm, covEm, covPm, estVarsZm )
	parsZfm   <- list( pathRg, CVfm, CVmf)
	modelMZf  <- mxModel( parsZf, meanGf, covMZf, dataMZf, expMZf, funML, name="MZf" )
	modelDZf  <- mxModel( parsZf, meanGf, covDZf, dataDZf, expDZf, funML, name="DZf" )
	modelMZm  <- mxModel( parsZm, meanGm, covMZm, dataMZm, expMZm, funML, name="MZm" )
	modelDZm  <- mxModel( parsZm, meanGm, covDZm, dataDZm, expDZm, funML, name="DZm" )
	modelDZo  <- mxModel( parsZf, parsZm, parsZfm, meanGfm, covDZo, dataDZo, expDZo, funML, name="DZo" )
	multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )
	QualAceModel  <- mxModel( "QualACE", parsZf, parsZm, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi )

	# ------------------------------------------------------------------------------
	# RUN MODEL

# Run Qualitative Sex Differences ACE model
QualAceFit    <- suppressWarnings(mxRun(QualAceModel, silent = T))
SUMM <- summary(QualAceFit)

Ests <-    as.data.frame(cbind(rbind(am, cm, af, cf, rg), 
	         rbind(QualAceFit$Am$result, QualAceFit$Cm$result, QualAceFit$Af$result, QualAceFit$Cf$result, QualAceFit$DZo$rg$values)))
colnames(Ests) <- c("Simulated", "Estimated")
Ests$Difference <- Ests$Simulated - Ests$Estimated
Imprecision <- sum(Ests$Difference)

# Fit submodels
CMfit <- omxSetParameters(QualAceFit, labels = "cm11", free = F, values = 0) ; CMfit <- suppressWarnings(mxRun(CMfit, silent = T))
CFfit <- omxSetParameters(QualAceFit, labels = "cf11", free = F, values = 0) ; CFfit <- suppressWarnings(mxRun(CFfit, silent = T))

RGfit <- omxSetParameters(QualAceFit, labels = "rg11", free = F, values = 1)   ; RGfit <- suppressWarnings(mxRun(RGfit, silent = T))

AMfit <- omxSetParameters(RGfit, labels = "am11", free = F, values = 0) ; AMfit <- suppressWarnings(mxRun(AMfit, silent = T) )
AFfit <- omxSetParameters(RGfit, labels = "af11", free = F, values = 0) ; AFfit <- suppressWarnings(mxRun(AFfit, silent = T))

EqualAfit  <- omxSetParameters(RGfit,      labels = c("am11", "af11"), newlabels = "a", free = T, values = 0) ; EqualAfit <- suppressWarnings(mxRun(EqualAfit, silent = T) )
EqualCfit  <- omxSetParameters(QualAceFit, labels = c("cm11", "cf11"), newlabels = "c", free = T, values = 0) ; EqualCfit <- suppressWarnings(mxRun(EqualCfit, silent = T) )
EqualEfit  <- omxSetParameters(QualAceFit, labels = c("em11", "ef11"), newlabels = "e", free = T, values = 0) ; EqualEfit <- suppressWarnings(mxRun(EqualEfit, silent = T) )

ncpAm <- (AMfit$output$Minus2LogLikelihood - RGfit$output$Minus2LogLikelihood)/(TotalN)
ncpAf <- (AFfit$output$Minus2LogLikelihood - RGfit$output$Minus2LogLikelihood)/(TotalN)

ncpCm <- (CMfit$output$Minus2LogLikelihood - QualAceFit$output$Minus2LogLikelihood)/(TotalN)
ncpCf <- (CFfit$output$Minus2LogLikelihood - QualAceFit$output$Minus2LogLikelihood)/(TotalN)

ncpRg <- (RGfit$output$Minus2LogLikelihood - QualAceFit$output$Minus2LogLikelihood)/(TotalN)


ncpEqualA <- (EqualAfit$output$Minus2LogLikelihood - RGfit$output$Minus2LogLikelihood)/(TotalN)
ncpEqualC <- (EqualCfit$output$Minus2LogLikelihood - QualAceFit$output$Minus2LogLikelihood)/(TotalN)
ncpEqualE <- (EqualEfit$output$Minus2LogLikelihood - QualAceFit$output$Minus2LogLikelihood)/(TotalN)


if( Imprecision > .02 ) warning("The estimated values may differ from the simulated values more than you are comfortable with. 
Make sure you double check these values. You may want to increase the sample size ", call. = F)
	
return(list(Parameters = Ests, WncpAm = ncpAm, WncpAf = ncpAf, WncpCm = ncpCm, WncpCf = ncpCf, WncpRg = ncpRg, 
	WncpEqualA = ncpEqualA, WncpEqualC = ncpEqualC, WncpEqualE = ncpEqualE))


	}


	# ------------------------------------------------------------------------------
	# Data Simulation Functions

	# Function to simulate Cholesky Twin Data
	CholSim <- function(aLower, cLower, eLower, Nmz, Ndz){
	
		nv <- dim(aLower)[1]
		ntv <- 2* nv

		II <- diag(2)
		UU <- matrix(1,2,2)
		GG <- matrix(c(1,.5,.5,1), 2, 2)

		impCovMZ <-   (UU %x% (aLower  %*% t(aLower)) + 
		               UU %x% (cLower  %*% t(cLower)) + 
		               II %x% (eLower  %*% t(eLower)) )  

		impCovDZ <-   (GG %x% (aLower  %*% t(aLower)) + 
		               UU %x% (cLower  %*% t(cLower)) + 
		               II %x% (eLower  %*% t(eLower)) )  

		mzData <- mvrnorm(Nmz, rep(0, ntv), impCovMZ, empirical = T)
		dzData <- mvrnorm(Ndz, rep(0, ntv), impCovDZ, empirical = T)

	colnames(mzData) <- colnames(dzData) <- paste(rep(paste(rep("t", nv), 1:nv, sep = "_"), 2), c(rep(1, nv), rep(2,nv)), sep = "_")
	
		return(list(mzData=mzData, dzData=dzData))
	}
	
	# Function to simulate Independent Pathway Twin Data
	indPathSim <- function(aLoad, cLoad, eLoad, specA, specC, specE, Nmz, Ndz){
	
		nv <- length(aLoad)
		ntv <- 2* nv
		loadA <- matrix(aLoad, nv,1)
		loadC <- matrix(cLoad, nv,1)
		loadE <- matrix(eLoad, nv,1)

		II <- diag(2)
		UU <- matrix(1,2,2)
		GG <- matrix(c(1,.5,.5,1), 2, 2)

		impCovMZ <-   (UU %x% (loadA  %*% t(loadA)) + 
		               UU %x% (loadC  %*% t(loadC)) + 
		               II %x% (loadE  %*% t(loadE)) +  
		               ( UU %x% (vec2diag(specA) %*% t(vec2diag(specA)))
		               + UU %x% (vec2diag(specC) %*% t(vec2diag(specC)))  
		               + II %x% (vec2diag(specE) %*% t(vec2diag(specE)))    )) 

		impCovDZ <-      (GG %x% (loadA  %*% t(loadA)) + 
		                  UU %x% (loadC  %*% t(loadC)) + 
		                  II %x% (loadE  %*% t(loadE)) +  
		                  ( GG %x% (vec2diag(specA) %*% t(vec2diag(specA)))
		                  + UU %x% (vec2diag(specC) %*% t(vec2diag(specC)))  
		                  + II %x% (vec2diag(specE) %*% t(vec2diag(specE)))    )) 

		mzData <- mvrnorm(Nmz, rep(0, ntv), impCovMZ, empirical = T)
		dzData <- mvrnorm(Ndz, rep(0, ntv), impCovDZ, empirical = T)
		colnames(mzData) <- colnames(dzData) <- paste(rep(paste(rep("t", nv), 1:nv, sep = "_"), 2), c(rep(1, nv), rep(2,nv)), sep = "_")
	
		return(list(mzData=mzData, dzData=dzData))
	}	
	
	# Function to simulate Common Pathway Twin Data
	comPathSim <- function(aL, cL, eL, load, specA, specC, specE, Nmz, Ndz){

	nv <- length(load)
	ntv <- 2* nv
	aLat <- matrix(sqrt(aL), 1,1)
	cLat <- matrix(sqrt(cL), 1,1)
	eLat <- matrix(sqrt(eL), 1,1)
	LOAD <- matrix(c(load), nv, 1)

	II <- diag(2)
	UU <- matrix(1,2,2)
	GG <- matrix(c(1,.5,.5,1), 2, 2)

	impCovMZ <- (II %x% LOAD) %*% ( UU %x% (aLat %*% t(aLat)) + UU %x% (cLat %*% t(cLat)) + II %x% (eLat %*% t(eLat))   )  %*% t(II %x% LOAD)+ 
	   ( UU %x% (vec2diag(specA) %*% t(vec2diag(specA)))
	   + UU %x% (vec2diag(specC) %*% t(vec2diag(specC)))  
	   + II %x% (vec2diag(specE) %*% t(vec2diag(specE)))    ) 

	impCovDZ <- (II %x% LOAD) %*% ( GG %x% (aLat %*% t(aLat)) + UU %x% (cLat %*% t(cLat)) + II %x% (eLat %*% t(eLat))   )  %*% t(II %x% LOAD)+ 
	   ( GG %x% (vec2diag(specA) %*% t(vec2diag(specA)))
	   + UU %x% (vec2diag(specC) %*% t(vec2diag(specC)))  
	   + II %x% (vec2diag(specE) %*% t(vec2diag(specE)))    ) 

	mzData <- mvrnorm(Nmz, rep(0, ntv), impCovMZ, empirical = T)
	dzData <- mvrnorm(Ndz, rep(0, ntv), impCovDZ, empirical = T)

	colnames(mzData) <- colnames(dzData) <- paste(rep(paste(rep("t", nv), 1:nv, sep = "_"), 2), c(rep(1, nv), rep(2,nv)), sep = "_")
	return(list(mzData=mzData, dzData=dzData))
	}





