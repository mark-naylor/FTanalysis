
getMinAgeMod_Verm <- function(nS, nI, doLog=TRUE) {
  
#   gspbest = array(dim=3)  # gsp = gamma, sigma, pi (->Age, dispersion, proportion at minimum)
  
  zUs = (nS+0.5)/(nI+0.5)
  if( doLog==TRUE){ zUs = log(zUs) }
  
  sUs = sqrt( (1/(nS+0.5)) + (1/(nI+0.5)) )
  nGrains = length(nI)
  
  print(zUs)
  print(sUs)
  
    ng = 100; ns = 25 ; np = 10    #  number of iterations for each parameter;
    
#     gUs = nS/nI ; gUs[which(gUs==0)]=0.0001 ; n=which(gUs==Inf) ; gUs[n]= (nS[n]+0.5)/(nI[n]+0.5)
    gUs=zUs
#     if(doLog) { gUs = log(gUs) } 
    ming = min(gUs)
    maxg = max(gUs)
    dg = (maxg-ming)/ng           # Increment for Gamma
    if(doLog) { ds = 1/ns } else { ds = maxg/ns }           # Increment for sigma
    dp = 1/np           # Increment for proportions
    gsp = c(ming,1,1)
    gspbest=gsp
    oldLL = LL(zUs, sUs, nGrains, gsp)
    # loop through gamma at higher resolution
    for (i in 1:ng){
      gsp[1] = ming + (i-1)*dg
      for (j in 1:ns){
        gsp[2] = j*ds
        for (k in 1:np){
          gsp[3] = k/(np+1)
          newLL = LL(zUs, sUs, nGrains, gsp)
          if (newLL>oldLL){
            gspbest <- gsp
            oldLL <- newLL
          }
        }
      }
    }
    print(paste( dg=dg/2, ds=ds/2, dp=dp/2, nGrains, sep="  ;  ") )
    gammaErr = cov_Verm(zUs, sUs, dg=dg/2, ds=ds/2, dp=dp/2, gspbest, nGrains)
  print(gammaErr)
print(gspbest)
#   } catch (Exception e){
#     gspbest[1] = tmM[1];
#     gammaErr = tmM[3];
#   }
  if (doLog){
    gspbest[2] = exp(gspbest[1]+gspbest[2])-exp(gspbest[1])  # ??? Ask Pieter - experr function...!
    gspbest[1] = exp(gspbest[1])
    gammaErr   = gammaErr*gspbest[1]
  }

print(gammaErr)
print(gspbest)
print(oldLL)

results = data.frame(gamma=gspbest[1] ,sigma=gspbest[2], prop=gspbest[3], gammaErr=gammaErr, logLike=oldLL)
  return(results)
}

LL <- function (zUs, sUs, nGrains, gsp){
  LL = 0
  for (u in 1:nGrains){
    zu = zUs[u]  ;  su = sUs[u]
    LL = LL + log( fu_Verm(zu, su, gsp) )
  }
  return( LL )
}

fu_Verm <- function (zu, su, gsp){
  gamma = gsp[1]
  sigma = gsp[2]
  prop = gsp[3]
  vu = su*su
  v = sigma*sigma
  
  mu0u = (gamma/v + zu/vu)/(1/v + 1/vu)
  s0u = 1/sqrt(1/v + 1/vu)
  A = prop/sqrt(2*pi*vu)
  B = -(zu-gamma)^2 /(2*vu)
  C = (1-prop) / sqrt( 2*pi*(v+vu) )
#   D = ( 1 - pnorm(mean=0, sd=1, q=(gamma-mu0u)/s0u)) / (1 - pnorm(mean=0, sd=1, q=0) )
  D = ( 1 - pnorm(mean=0, sd=1, q=(gamma-mu0u)/s0u)) * 2
  E = -(zu-gamma)^2 / (2*(v+vu))
  return( A*exp(B) + C*D*exp(E) )
}

cov_Verm <- function (zUs, sUs, dg, ds, dp, gspbest, nGrains) {
#   dg=dg/2; ds=ds/2; dp=dp/2
  n = nGrains
  g = gspbest[1]  ;  s = gspbest[2]  ;  p = gspbest[3]
  J = array(dim=c(3,3)) # Jacobian
  J[1,1] = -(LL(zUs, sUs,n,c(g+dg,s,p))-2*LL(zUs, sUs,n,gspbest)+LL(zUs, sUs,n,c(g-dg,s,p)))/(dg*dg);
  J[2,2] = -(LL(zUs, sUs,n,c(g,s+ds,p))-2*LL(zUs, sUs,n,gspbest)+LL(zUs, sUs,n,c(g,s-ds,p)))/(ds*ds);
  J[3,3] = -(LL(zUs, sUs,n,c(g,s,p+dp))-2*LL(zUs, sUs,n,gspbest)+LL(zUs, sUs,n,c(g,s,p-dp)))/(dp*dp);
  J[1,2] = -(LL(zUs, sUs,n,c(g+dg,s+ds,p))-LL(zUs, sUs,n,c(g+dg,s-ds,p))-LL(zUs, sUs,n,c(g-dg,s+ds,p))+LL(zUs, sUs,n,c(g-dg,s-ds,p)))/(4*dg*ds);
  J[1,3] = -(LL(zUs, sUs,n,c(g+dg,s,p+dp))-LL(zUs, sUs,n,c(g+dg,s,p-dp))-LL(zUs, sUs,n,c(g-dg,s,p+dp))+LL(zUs, sUs,n,c(g-dg,s,p-dp)))/(4*dg*dp);
  J[2,3] = -(LL(zUs, sUs,n,c(g,s+ds,p+dp))-LL(zUs, sUs,n,c(g,s+ds,p-dp))-LL(zUs, sUs,n,c(g,s-ds,p+dp))+LL(zUs, sUs,n,c(g,s-ds,p-dp)))/(4*ds*dp);
  J[2,1] = J[1,2];
  J[3,1] = J[1,3];
  J[3,2] = J[2,3];
  det = J[1,1]*J[2,2]*J[3,3] + J[2,1]*J[3,2]*J[1,3] + J[3,1]*J[1,2]*J[2,3] -
    J[1,1]*J[3,2]*J[2,3] - J[3,1]*J[2,2]*J[1,3] - J[2,1]*J[1,2]*J[3,3];
  gammaVar = (J[3,3]*J[2,2] - J[2,3]*J[3,2])/det;
  try({  gammaErr = sqrt(gammaVar)  })
  
  return(gammaErr)
}

