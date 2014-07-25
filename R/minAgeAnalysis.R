library(pracma)

fU4 <- function(zU, sigmaU, gamma, mu, sigma, piGamma){
  mu0U = (mu/sigma**2  +  zU/sigmaU**2) / (1/sigma**2  +  1/sigmaU**2)
  sigma0U = 1/ sqrt( (1/sigma**2) + (1/sigmaU**2) )
  
  p = (gamma-mu)/sigma/sqrt(2)
  p0 = (gamma-mu0U)/sigma0U/sqrt(2)
  Phi = 0.5 * (1+ erf(p))
  Phi0 = 0.5 * (1+erf(p0))
  
  fU = piGamma/sqrt(2*pi*sigmaU**2) * exp(- 0.5 * (zU-gamma)**2 / sigmaU**2)  +
    (1-piGamma) /sqrt(2*pi*(sigmaU**2+sigma**2)) * (1-Phi0) / (1-Phi) * exp(- 0.5 * (zU-mu)**2 / (sigmaU**2+sigma**2)) 
  
  Re( sum( log(fU) ) )
}


fU3 <- function(zU, sigmaU, gamma, sigma, piGamma){
  mu<-gamma
  mu0U = (mu/sigma**2  +  zU/sigmaU**2) / (1/sigma**2  +  1/sigmaU**2)
  sigma0U = 1/ sqrt((1/sigma**2) + (1/sigmaU**2))
  
  p0 = (gamma-mu0U)/sigma0U/sqrt(2)
#    if(p0 > 10){  
#      Phi0 = 1 
#    } else if(p0 < -10){  
#      Phi0 = 0
#    } else { 
    Phi0 = 0.5 * (1+erf(p0)) 
#    }

  
  fU = piGamma/sqrt(2*pi*sigmaU**2) * exp(- 0.5 * (zU-gamma)**2 / sigmaU**2)  +
    (1-piGamma) /sqrt(2*pi*(sigmaU**2+sigma**2)) * (1-Phi0) / (0.5) * exp(- 0.5 * (zU-mu)**2 / (sigmaU**2+sigma**2)) 
  
  Re( sum( log(fU) ) )
}


findMinAgeGaussian3 <- function(piGammaRange=seq(0.1,1,0.1), sigmaRange=seq(0.1,5,0.1), ageRange, FTdataset){
  lambdaD = 1.55125E-4
  nI=FTdataset$nI
  nS=FTdataset$nS
  zeta=FTdataset$Zeta
  rhoD=FTdataset$rhoD/1E6
  
  zU = log((nS+0.5)/(nI+0.5))
  sigmaU = sqrt( (1/(nS+0.5)) + (1/(nI+0.5)) )
  
  trackRatioRange = (exp(ageRange*lambdaD)-1)/(0.5*lambdaD*zeta*rhoD)
  gammaRange = log(trackRatioRange)
  
  gridExtent3 <- expand.grid(
    piGamma=piGammaRange,
    gamma=gammaRange,
    sigma=sigmaRange
  )
  
  L<-c()
  for(i in 1:length(gridExtent3[,1])) {
    L[i] <- fU3(zU, sigmaU, gridExtent3$gamma[i], gridExtent3$sigma[i], gridExtent3$piGamma[i]) 
  }
  
  Lmax=max(L, na.rm=TRUE)
  n=which(L==Lmax )
  results = gridExtent3[n,]
  
  trackRatio = exp(results$gamma)
  
  age = 1/lambdaD * log(1+0.5*lambdaD*zeta*rhoD* trackRatio)
  
  c(age, gamma, sigma, piGamma, L[n])
}
