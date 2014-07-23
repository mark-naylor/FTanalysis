library(NORMT3)

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
  Phi0 = 0.5 * (1+erf(p0))
  
  fU = piGamma/sqrt(2*pi*sigmaU**2) * exp(- 0.5 * (zU-gamma)**2 / sigmaU**2)  +
    (1-piGamma) /sqrt(2*pi*(sigmaU**2+sigma**2)) * (1-Phi0) / (0.5) * exp(- 0.5 * (zU-mu)**2 / (sigmaU**2+sigma**2)) 
  
  Re( sum( log(fU) ) )
}
