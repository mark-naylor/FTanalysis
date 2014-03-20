# '... BinomFit.BAS  v. 1.8  9/15/2004   Mark Brandon, Yale University
# '============================================================================
# '... Determines best-fit peaks based on an initial guess by the user
# ' using the method of Galbraith and Green (1990, Nuclear Tracks and
# ' Radiation Measurements, v. 17, p. 197-206). The best-fit peaks are
# ' found using the EM method as outlined in Galbraith (1988, Technometrics,
# ' v. 30, p. 280).
# '... 2/19/92  The binomial probability density algorithm halted with
# ' an overflow error when the total number of tracks was large. Thus, Function
# ' FIU was calculated using a Gaussian approximation when PKtheta was between
# ' 5% and 95% and Ru% > 100. This modification increase the speed of the
# ' program significantly. The PKtheta bounds are need to adhere to the
# ' Gaussian approximation. If overflow continues to occur then it might be
# ' necessary to introduce a Poisson approximation as well.
# '... 6/7/92 version 1.3  Standard error for peak age is now calculated using
# ' the method of Galbraith (Nuclear Tracks, 1990, eqns. 6, 7), which produces
# ' asymmetric and more accurate confidence limits.
# '... 2/1/93 version 1.4  Minor changes so that grain ages and initial guess
# ' for the probability mass of the peaks are calculated using Z method of
# ' Galbraith and Green (1990, Nuclear Tracks, v. 17, p. 197-206.)
# '... 11/10/93 version 1.5  Includes a test for non-convergence of a solution
# ' in the EM subroutine.
# '... 11/28/94 version 1.6  Now reports SD(z), which is approximately equal
# ' to W, the peak width in a PD plot. SD(z) is estimated assuming a Gaussian
# ' distribution for z with the mean and std dev estimated from the Binomial
# ' distribution of Theta.  Also introduced an improved version DiskOpen.
# '... 3/12/95 version 1.7  Fixed minor bug associated with Z method and
# ' modified the input routine to ignore blank and comment lines.
# '... 6/24/97 version 1.8  Minor change to calculation of estimated
# ' peak width, and extended reporting of uncertainties.
# '... 2/5/98: Included additional note in introduction of program.
# '... 4/18/01: Fixed a rare division-by-zero problem in the EM routine.
# '... 6/11/01: Fixed header for 68% confidence interval; note that no
# ' caluclation problems were involved.
# '... 3/04: Added informal option at the end of EM subroutine to write probabilities
# ' for grain ages in a separate file.
# '... 9/15/04: The division-by-zero problem was not completely solved with
# ' the previous fix (4/18/01). I made a more complete fix by imposing a lower
# ' probability limit for Fiu and a lower limit for PkFrac(I%)for 0.001 in 
# ' PkFind in order to avoid round-off errors.
# MN
# Remove FactLn in favour of lfactorial()
# GammLn Then redundant
# '============================================================================

library(zipfR)

LamdaD <- 1.55125E-10  # 'Total decay constant for 238U (yr^-1)
GF <- .5  
C0 <- .3989423#        # 'Used in Gaussian equation, equals 1/sqrt(2*Pi)

Bico <- function (N, K){
  # '============================================================================
  # '... Returns the binomial coefficient (n,k) as a floating-point number.
  # ' From Press and others (1986, p. 157).
  # '============================================================================
  d = exp( (lfactorial(N) - lfactorial(K) - lfactorial(N - K)) )
  Bico = d
  Bico
}


Chi2 <- function (PkTheta, PkFrac, PkNum, Num, Nt, Ns) {
  # '============================================================================
  # '... Calculates Chi^2 parameter for the binomial model of Galbraith
  # ' and Green. PkTheta() contains the theta ratios for the best-fit peaks, and
  # ' Pkfrac() contains the proportions of the total distribution for the peaks.
  # ' Num% (shared) = number of dated grains
  # ' PkNum% (shared) = number of peaks to be found
  # ' Ns() (shared) = number of spontaneous tracks counted per grain
  # ' Nt() (shared) = total number of tracks counted per grain
  # ' Chi2 = function is returned with Chi^2 value
  # '============================================================================
  # '... Find mean and sum of squares for solution vector
  Mu = 0
  SS = 0
  for (I in 1:PkNum){
    Mu = Mu + PkFrac[I] * PkTheta[I]
    SS = SS + PkFrac[I] * PkTheta[I] * PkTheta[I]
  }
  
  # '... Find Chi2 misfit for this particular solution vector
  Factor1 = Mu * (1 - Mu)
  Factor2 = SS - Mu * Mu
  C2 = 0
  for (I in 1:Num){
    NNt = Nt[I]
    dev = Ns[I] - Mu * NNt
    Var = Factor1 * NNt + NNt * (NNt - 1) * Factor2
    C2 = C2 + dev * dev / Var
  }
  Chi2 = C2
  return(Chi2)
}

EM <- function(PkTheta, PkFrac, PkNum, Num, Ns, Ni, Nt, verbose=TRUE){
  # '============================================================================
  # '... Uses the EM method outlined in Galbraith (1988) to find the best-fit
  # ' peaks (see Titterington and others, 1985, p. 84-97 for further details).
  # ' Convergence is defined as when the relative changes in the Log-Likelihood
  # ' value are less than the constant Rtol.
  # '============================================================================
  Rtol = .00001
  P=array(dim=c(10, 230))
  # '============================================================================
  
  if(verbose){ print(  "ITERATION  LOG-LIKELIHOOD  RELATIVE CHANGE" ) }
  
  OldLogLike = 0
  for(Iter in 1:300){
    # '... Calculate new values for the array P() and LogLikelihood parameter
    LogLike = 0
    for (U in 1:Num) {
      SumFiu = 0
      for (I in 1:PkNum){
        P[I, U] = PkFrac[I] * Fiu(Ns[U], Ni[U], Nt[U], PkTheta[I])
        SumFiu = SumFiu + P[I, U]
        # if (P[I, U] == 0) { INPUT Q }
      }
      for (I in 1:PkNum){
        P[I, U] = P[I, U] / SumFiu
      }
      LogLike = LogLike + log(SumFiu)
    }
    
    # '... Update proportions of the peaks in Pkfrac()
    for (I in 1:PkNum){
      SumPkfrac = 0
      for (U in 1:Num){
        SumPkfrac = SumPkfrac + P[I, U]
      }
      PkFrac[I] = SumPkfrac / Num
    }
    # '... Update peak theta ratios in PkTheta()
    for (I in 1:PkNum){
      Sum1 = 0
      Sum2 = 0
      for (U in 1:Num){
        Sum1 = Sum1 + P[I, U] * Ns[U]
        Sum2 = Sum2 + P[I, U] * Nt[U]
      }
      PkTheta[I] = Sum1 / Sum2
    }
    
    ChiSq = Chi2(PkTheta, PkFrac, PkNum, Num, Nt, Ns)
    # '... Check for convergence: the relative size of the final steps must
    # ' remain small, less than Rtol. The routine now includes a limit
    # ' on the maximum number of iterations that can occur without advancing
    # ' the solution.
    Rchange = abs(.5 * (LogLike - OldLogLike) / (LogLike + OldLogLike))
    if(verbose){ print( paste( Iter, LogLike, Rchange ) ) }
    
    if (Rchange <= Rtol) {
      Itest1 = Itest1 + 1
      if (Itest1 == 3) { break }
    } else {
      Itest1 = 0
    }
    OldLogLike = LogLike
    
    if (Iter > 300){
      print(  "   >>>    WARNING: LOCAL MINIMUM POORLY DEFINED   <<<" )
      print(  "   >>> CHECK THAT BEST-FIT SOLUTION IS REASONABLE <<<" )
      break
    }
  }
  
  
  
  # OPEN "Binomfit.prb" for OUTPUT AS #9
  # for (U in 0:Num){
  # if (U == 0) {
  # print(   Label1 )
  # print(   Label2 )
  # print(   "FILE HEADER:" )
  # print(   Label3 )
  # print(  "" )
  # print(   "Probability in precent for best fit, on a grain and peak basis" )
  # }
  # for (I in 0:PkNum){
  
  # if (U == 0){
  # if (I == 0) {
  # print(  #9, "Grain# ";)
  # } else {
  # print(  #9, USING "\      \###   "; "Prob_Pk#"; I%;)
  # }
  # } else {
  # if (I == 0) {
  # print(  #9, USING "####   "; U; )
  # } else {
  # print(  #9, USING "     ###.###  "; P[I, U] * 100;)
  # }
  # }
  # }
  # print(  #9, )
  # }
  
  
  results = data.frame( PkTheta, PkFrac, LogLike, Iter)
  return(results)
}



Fiu <- function (Ru, Su, Tu, PkTheta){
  # '============================================================================
  # '... Evaluates the contribution that the Uth grain makes to the binomial
  # ' probability density for the Ith peak. The binomial probability density
  # ' calculation halts with an overflow error for Ru% greater than about 500.
  # ' As a result, the Gaussian probability density distribution is used to
  # ' approximate the binomial probability density function for Ru% > 100 and
  # ' PKtheta between 5% and 95% (see Taylor, 1982, p. 194, for details).
  # ' Ru%, Su%, Tu% = spontaneous, induced, and total track for the Uth grain.
  # ' PkTheta = mean theta ratio for the Ith peak.
  # ' The constant C0 = 1/SQR(2*Pi) and is set in the main program.
  # ' ... 9/15/04: added option to ensure a minimum value for Fiu, in order
  # ' to avoid problems with round-off errors.
  # '============================================================================
  # # # print(paste("a",Tu, "b",PkTheta, "c",Ru, "d",Su))
  if ((Tu > 100) && (PkTheta > .05) && (PkTheta < .95)){
    Mean = Tu * PkTheta
    SD = sqrt(Mean * (1 - PkTheta))
    Z = (Ru - Mean) / SD
    Tmp = (C0 / SD) * exp(-0.5 * Z * Z)
  } else {  
    Tmp = Bico(Tu, Ru) * (PkTheta ^ Ru) * (1 - PkTheta) ^ Su
  }
  
  if (Tmp <= 1E-10){ Fiu = 1E-10 } else {Fiu = Tmp}
  
  Fiu
}



# # PkFind <-function(PeakAge, PkTheta, PkFrac){
# # '============================================================================
# # '... Asks the user for an initial guess of the peaks in the observed
# # ' composite probability density plot.
# # ' PeakAge = age of selected peaks
# # ' PkTheta = theta ratio for each peak
# # ' Pkfrac = fractional proportion of the total distribution of each peak
# # '        where Sum[Pkfrac()] = 1
# # ' Num% (shared) = number of grain ages or peaks in the input dataset
# # ' PkNum% (shared) = total number of peaks
# # ' Zgrain(), ZErr() (shared) = Z and SE(Z) for grains
# # ' AgeMin!, AgeMax! (shared) = age range for grain-age distribution
# # '============================================================================
# # '... Set up peak list
# print (paste("Number of peaks to be fitted: ", PkNum ))

# for (I in 1:PkNum){
# PeakAge[I] = runif(1, AgeMin, AgeMax)
# PkTheta[I] = Theta(PeakAge[I])
# }

# # '... Make a guess of proportional size of each peak
# for (I in 1:PkNum){
# # '... Calculate probability density at the peak.
# Zpk = ZfromTau( PeakAge[I] )
# PkDen = SumPD(Zpk, 0)
# # '... Fraction of total grains = peak density * std dev *sqrt(2*Pi)/ Num%
# # ' Use a SD(z) = 0.15 to estimate the proportion of each peak
# # ' relative to the total distribution.
# PkFrac[I] = PkDen * .15 / (C0 * Num)
# # '... Set a lower limit for PkFrac(I%) to avoid division-by-zero problems.
# if (PkFrac[I] < .001) { PkFrac[I] = .001 }
# }
# }

StdErr <- function (FTdataset , PkTheta, PkFrac, PkNum, GF){
  # FTdataset=FTdataset; PkTheta=results$PkTheta; PkFrac=results$PkFrac; PkNum=PkNum; GF=GF	
  # '============================================================================
  # '... Calculates peak ages and approximate standard errors for the best-fit
  # ' parameters using the method outlined in Galbraith (1988, p. 280). Note that
  # ' the method he uses produces a scaled covariance matrix.
  # ' ZErr() (shared) = SE(Z) for grains
  # '============================================================================
  # '... Calculate the array P()
  
  Ns <- FTdataset$nS
  Ni <- FTdataset$nI
  Nt <- FTdataset$nT
  Num <- FTdataset$nGrain
  Zeta0 <- FTdataset$Zeta
  REZeta0  <- FTdataset$relErrZeta
  RhoD0 <- FTdataset$rhoD
  RERhoD0 <- FTdataset$relErrRhoD
  
  P=array(dim=c(10, 230))
  A=array(dim=c(9, 9))
  B=array(dim=c(9, 10))
  C=array(dim=c(10, 10))
  Covar=array(dim=c(19, 19))
  SEPkfrac = array( dim=(10))
  PeakAge = array( dim=( 10 ) )
  SEPeakAge = array( dim=c(10, 4))
  PkSDz = array( dim=(10) )
  
  for (U in 1:Num){
    SumFiu = 0
    for (I in 1:PkNum){
      P[I, U] = PkFrac[I] * Fiu(Ns[U], Ni[U], Nt[U], PkTheta[I])
      SumFiu = SumFiu + P[I, U]
    }
    for (I in 1:PkNum){
      P[I, U] = P[I, U] / SumFiu
    }
  }
  
  
  # '... Construct A(), the upper right quadrant of the Inverse Covariance matrix
  PIk = PkFrac[PkNum]
  for (I in 1:(PkNum-1)){
    for ( J in 1:(PkNum-1)){
      A[I, J] = 0
      for (U in 1:Num){
        Pku = P[PkNum, U]
        Term = (P[I, U] / PkFrac[I] - Pku / PIk) * (P[J, U] / PkFrac[J] - Pku / PIk)
        A[I, J] = A[I, J] + Term
      }
    }
  }
  # '... Construct B(), the diagonal quadrants of the Inverse Covariance matrix
  for (I in 1:(PkNum-1)){
    for (J in 1:PkNum){
      B[I, J] = 0
      for (U in 1:Num){
        Pku = P[PkNum, U]
        if(I == J){  Del1 = 1 } else { Del1 = 0 }
        if(PkNum == J) { Del2 = 1 } else { Del2 = 0 }
        Term = P[I, U] / PkFrac[I] - Pku / PIk - Del1 / PkFrac[I] + Del2 / PIk
        Aju = Ns[U] - PkTheta[J] * Nt[U]
        B[I, J] = B[I, J] + P[J, U] * Aju * Term
      }
    }
  }
  # '... Construct C(), the lower left quadrant of the Inverse Covariance matrix
  for (I in 1:PkNum){
    for (J in 1:PkNum){
      C[I, J] = 0
      for (U in 1:Num){
        # '... Calculate first derivatives of log(fiu) wrt beta(i)
        # ' and log(fju) wrt beta(j)
        Aiu = Ns[U] - PkTheta[I] * Nt[U]
        Aju = Ns[U] - PkTheta[J] * Nt[U]
        Term = P[J, U] * Aiu * Aju
        if (I == J){
          # '... Add second derivative of log(fiu) wrt beta(i)
          Term = Term - (Aiu * Aiu - PkTheta[I] * (1 - PkTheta[I]) * Nt[U])
        }
        C[I, J] = C[I, J] + P[I, U] * Term
      }
    }
  }
  
  # '... Assembly inverse of covariance matrix
  for (I in 1:(PkNum-1)){
    for (J in 1:(PkNum-1)){
      Covar[I, J] = A[I, J]
    }
  }
  for (I in 1:(PkNum-1)){
    for (J in 1:PkNum){
      Covar[I, J + PkNum - 1] = B[I, J]
    }
  }
  for (I in 1:PkNum){
    for (J in 1:(PkNum-1)){
      Covar[I + PkNum - 1, J] = B[J, I]
    }
  }
  for (I in 1:PkNum){
    for (J in 1:PkNum){
      Covar[I + PkNum - 1, J + PkNum - 1] = C[I, J]
    }
  }
  
  # '... Invert the matrix Covar using the SVD routine
  M = 2 * PkNum - 1
  
  # for(i in 1:M){
  # for(j in 1:M){
  # Covar[i,j] <- Covar[i,j] + (runif(1)-0.5)*2*0.000001
  # }
  # }
  # print(Covar)
  
  Covar[1:M,1:M] <- solve(Covar[1:M,1:M])
  # CALL SvdCmp(Covar(), W(), V(), M%, M%)
  # CALL SvdInv(Covar(), W(), V(), Cnum, M%)
  # print(  "Condition number for inversion of covariance matrix: "; Cnum
  
  # '... Calculate standard errors for the fraction parameter for 1 to PkNum%-1
  for (I in 1:(PkNum-1)){
    SEPkfrac[I] = sqrt( abs(Covar[I, I]) )
  }
  
  # '... Calculate standard error for the last fraction parameter (see Press
  # # # # # # # # # # # # # # # # # # # # # ' et al., 1992, p. 692, for details about this calculation) 
  if(PkNum>1){
    Sum = 0
    for (I in 1:(PkNum-1)){
      for (J in 1:(PkNum-1)){
        print(I,J)
        Sum = Sum + Covar[I, J]
      }
    }
    SEPkfrac[PkNum] = sqrt(abs(Sum))
  }
  
  # '... Calculate peak ages, standard errors, and estimated peak widths
  for (I in 1:PkNum){
    BetaPk = log(PkTheta[I] / (1 - PkTheta[I]))
    PeakAge[I] = TaufromBeta(BetaPk, Zeta0, GF, RhoD0)
    J = I + PkNum - 1
    SE = sqrt( abs(Covar[J, J]) + RERhoD0 * RERhoD0 + REZeta0 * REZeta0)
    SEPeakAge[I, 1] = TaufromBeta(BetaPk - SE, Zeta0, GF, RhoD0) - PeakAge[I]
    SEPeakAge[I, 2] = TaufromBeta(BetaPk + SE, Zeta0, GF, RhoD0) - PeakAge[I]
    SEPeakAge[I, 3] = TaufromBeta(BetaPk - 1.96 * SE, Zeta0, GF, RhoD0) - PeakAge[I]
    SEPeakAge[I, 4] = TaufromBeta(BetaPk + 1.96 * SE, Zeta0, GF, RhoD0) - PeakAge[I]
    SumNt = 0
    SumFiu = 0
    for (U in 1:Num){
      SumNt = SumNt + P[I, U] * Nt[U]
      SumFiu = SumFiu + P[I, U]
    }
    SumNt = SumNt / SumFiu
    PkSDz[I] = sqrt(1 / (SumNt * PkTheta[I] * (1 - PkTheta[I])))
  }
  
  # '... Reclaim temporary array space
  # ERASE P, A, B, C, Covar, W, V
  
  Cnum <- -999
  
  PeakAgeCI65min=SEPeakAge[, 1]
  PeakAgeCI65plus=SEPeakAge[, 2]
  PeakAgeCI95min=SEPeakAge[, 3]
  PeakAgeCI95plus=SEPeakAge[, 4]
  results = data.frame(PeakAge, PeakAgeCI65min, PeakAgeCI65plus, PeakAgeCI95min, PeakAgeCI95plus, SEPkfrac, PkSDz, Cnum)
  
  return(results)
  
}



SumPD <- function (Zi, Order, Zgrain, Zerr, Num, KernelFactor = .6 ){
  # '============================================================================
  # '... Calculates the probability density function and its derivatives
  # ' at the specified value of Zi.
  # ' Shared variables: Zgrain, Zerr, Num
  # '============================================================================
  
  Dsum = 0
  for (I in 1:Num){
    dev = Zi - Zgrain[I]
    S = KernelFactor * Zerr[I]
    if (abs(dev) < 5 * S){
      Dev2 = dev * dev
      S2 = S * S
      ExpTerm = exp(-Dev2 / (2 * S2))
      if(Order==0)	{
        Dsum = Dsum + (C0 / S) * ExpTerm
      } else if (Order==1){
        S3 = S * S2
        Dsum = Dsum - (C0 * dev / S3) * ExpTerm
      } else if (Order==2){
        S5 = S2 * S2 * S
        Dsum = Dsum + (C0 * (Dev2 - S2) / S5) * ExpTerm
      } else if (Order==3){
        S7 = S2 * S2 * S2 * S
        Dsum = Dsum + (C0 * dev * (-Dev2 + 3 * S2) / S7) * ExpTerm
      }
    }
  }
  SumPD = Dsum
  SumPD
}


TaufromBeta <- function (Beta, Zeta0, GF, RhoD0) {
  # '============================================================================
  # '... This routine calculates the FT age from Beta = ln(Ns/Ni).
  # ' TaufromBeta = FT age in Ma
  # ' Zeta0 (shared) = zeta factor
  # ' RhoD0 (shared) = monitor track density used for this dataset
  # ' LamdaD (shared) = total decay constant for 238U (yr^-1)
  # ' GF (shared) = geometry factor
  # '============================================================================
  if (Beta > 80) { 
    BBeta <- 80 
  } else if (Beta < -80){
    BBeta <- -80
  } else {
    BBeta <- Beta
  }
  ExpZm = LamdaD * Zeta0 * GF * RhoD0 * exp(BBeta)
  TaufromBeta = log(1. + ExpZm) / (1000000. * LamdaD)
  TaufromBeta
}


TaufromZ <- function(Z) {
  # '============================================================================
  # '... Converts to age in Myr from Z.
  # ' LamdaD is a constant set in the main program.
  # '============================================================================
  TaufromZ = log(1 + exp(Z)) / (1000000 * LamdaD)
  TaufromZ
}


Theta <- function (age, Zeta0, GF, RhoD0){
  # '============================================================================
  # '... Given a FT age, this routine calculates Theta = Ns /(Ns+Ni).
  # ' Age = FT age in Ma
  # ' Zeta0 (shared) = zeta factor
  # ' RhoD0 (shared) = monitor track density used for this dataset
  # ' LamdaD (shared) = total decay constant for 238U (yr^-1)
  # ' GF (shared) = geometry factor
  # '============================================================================
  ExpTL = exp( age * 1000000 * LamdaD )
  Theta = (ExpTL - 1) / (LamdaD * Zeta0 * GF * RhoD0 + ExpTL - 1)
  Theta
}

ZfromTau <- function (Tau){
  # '============================================================================
  # '... Converts to Z from age in Myr.
  # ' LamdaD is a constant set in the main program.
  # '============================================================================
  ZfromTau = log(exp(Tau * 1000000 * LamdaD) - 1)
  ZfromTau
}


# '... Chi2Comp  v 1.1  11/12/93   Mark Brandon, Yale University
# '===========================================================================
# ' Compares X^2 values to determine if fit is significantly
# ' improved. See Bevington (1969, p. 200) for description of method.
# ' The routine for calculating the F distribution probability is from Press
# ' and others (1986).
# '===========================================================================
# DECLARE FUNCTION Betai! (A!, B!, X!)
# DECLARE FUNCTION Betacf! (A!, B!, X!)
# DECLARE FUNCTION GammLn! (XX!)
# CONST True = -1, False = 0

getChi2Comp <- function(Chi1, DF1 , Chi2, DF2){
  
  print( "===================Chi2Comp Program v 1.1  (Brandon 11/12/93)===================")
  print(  "Compares X^2 values to determine if fit is significantly")
  print(  "improved. See Bevington (1969, p. 200) for description of method.")
  print(  "The routine for calculating the F distribution probability is from Press")
  print(  "and others (1986).")
  print(  "  ")
  print(  "F = [Chi^2(1) - Chi^2(2)] / [Chi^2(2)/DF2]")
  print(  "    where 1 refers to the first fit and 2,")
  print(  "    the second fit with Chi^2(1) >= Chi^2(2).")
  print( "  ")
  print(  "The 'second fit becomes first fit' option allows the user to efficiently")
  print(  "analyze a series of results where X^2 is successively reduced due to")
  print(  "an increase in the number of fit parameters. This option is particularly")
  print(  "useful for examining data from the BINOMFIT AND GAUSSFIT programs.")
  print(  "To halt program, use Ctrl-Break." )
  
  print( paste ( "First: Chi^2, degrees of freedom: ", Chi1, DF1 ) )
  print( paste ( "Second: Chi^2, degrees of freedom: ", Chi2, DF2 ) )
  
  F = (Chi1 - Chi2) / (Chi2 / DF2)
  
  print ( paste( "F = ", F) )
  DF1 = DF1 - DF2
  X = DF2 / (DF2 + DF1 * F)
  
  Prob = Ibeta(a=DF2 / 2, b=DF1 / 2, x=X)
  
  print( paste( "Probability that F is by chance alone (%) ", 100 * Prob ))
  
  return( list(F=F, P=Prob) )
}


KernelPD <- function (Zi, Order, Zgrain, Zerr, Num, KernelFactor = .6){
  # '============================================================================
  # '... Calculates the probability density function and its derivatives
  # ' at the specified value of Zi. When Order = -1, the routine returns the
  # ' estimated probability density and standard error. Otherwise SEPD == 0.
  # ' Constant: KernelFactor, which is set in main program
  # '============================================================================
  PD = 0
  SSPD = 0
  for (I in 1:Num){
    Dev = Zi - Zgrain[I]
    s = KernelFactor * Zerr[I]
    if (abs(Dev) < (5 * s)) {
      Dev2 = Dev * Dev
      S2 = s * s
      ExpTerm = exp(-Dev2 / (2 * S2))
      if(Order == -1){
        Term = (C0 / s) * ExpTerm
        PD = PD + Term
        SSPD = SSPD + Term * Term
      } else if(Order == 0){
        PD = PD + (C0 / s) * ExpTerm
      } else if(Order == 1){
        S3 = s * S2
        PD = PD - (C0 * Dev / S3) * ExpTerm
      } else if(Order == 2){
        S5 = S2 * S2 * s
        PD = PD + (C0 * (Dev2 - S2) / S5) * ExpTerm
      } else if(Order == 3){
        S7 = S2 * S2 * S2 * s
        PD = PD + (C0 * Dev * (-Dev2 + 3 * S2) / S7) * ExpTerm
      }
    }
  }
  if ((Order == -1) && (Num > 1)){
    # '... SE of estimated probability density; formula is after
    # ' eqn. 3.7 in Silverman (1986, p. 36). Note that the formula
    # ' accounts for the fact that PD has units of number of grains
    # ' per Num% grains per unit Z.
    SEPD = sqrt( (Num * SSPD - PD * PD) / ( Num - 1))
  } else {
    SEPD = 0
  }
  results <- c(PD, SEPD)
  results
}


calcPooledAge <- function(nS, nI, Zeta, ZetaStErr, RhoD, nD){
  trackRatioSum = sum(nS)/sum(nI)
  pooledAge = (1/LamdaD) * log(LamdaD*trackRatioSum*RhoD*Zeta*0.5 +1)/1E6
  pooledAgeError = pooledAge * sqrt( 1/sum(nI) + 1/sum(nS) + 1/nD + (ZetaStErr/Zeta)**2)
  
  output = list(pooledAge= pooledAge, StErrPooledAge= pooledAgeError)
  return(output)
}


# "===================Chi2Comp Program v 1.1  (Brandon 11/12/93)==================="
# "Compares X^2 values to determine if fit is significantly"
# "improved. See Bevington (1969, p. 200) for description of method."
# "The routine for calculating the F distribution probability is from Press"
# "and others (1986)."
#
# "F = [Chi^2(1) - Chi^2(2)] / [Chi^2(2)/DF2]"
# "    where 1 refers to the first fit and 2,"
# "    the second fit with Chi^2(1) >= Chi^2(2)."


Ftest <- function(Chi1, Chi2, DF1, DF2){
  if(Chi2>Chi1){ tmp1 <- Chi1 ; tmp2 <- DF1 ; Chi1 <- Chi2 ; DF1 <- DF2 ; Chi2 <- tmp1 ; DF2<- tmp2 }
  
  F <- (Chi1 - Chi2) / (Chi2 / DF2)
  DFdiff <- DF1 - DF2
  X <- DF2 / (DF2 + DFdiff * F)
  Prob <- Ibeta(a=DF2 / 2, b=DFdiff / 2, x=X)
  return( list(F=F, P=Prob) )
}



getPDFunction <- function(FTdataset, zeroNsOffset=0) {
  
  if(zeroNsOffset==0){
    Zgrain = FTdataset$Zgrain
    Zerr = FTdataset$Zerr
  } else {
    b <- FTdataset$b
    NNs <- FTdataset$nS + zeroNsOffset
    NNi <- FTdataset$nI + zeroNsOffset
    Zgrain = log(b * NNs / NNi)
    Zerr = sqrt(1 / NNs + 1 / NNi)	
  }
  
  Zmax = max(Zgrain)
  Zmin = min(Zgrain)
  
  Num = length(Zgrain)
  QMeanZerr = sum( Zerr * Zerr )
  QMeanZerr = sqrt(QMeanZerr / (Num))
  
  barWidth=0.1			# 'width of histogram bar in Z units
  zWidth = barWidth / 5	# 'interval width for PD plot in Z units
  
  LwLmt = barWidth * (floor(Zmin / barWidth) - 1)
  UpLmt = barWidth * (1 + floor(Zmax / barWidth))
  
  # '... Calculate probability density distribution and write output data.
  # '... Start and end at 5*QMeanZerr below Zmin and above Zmax, respectively.
  LwLmt = Zmin - (5 * QMeanZerr)
  UpLmt = Zmax + (5 * QMeanZerr)
  Zi = LwLmt
  atZ=c()
  # '  set above: zWidth = BarWidth / 5
  gCount = floor((UpLmt - LwLmt) / zWidth) + 1
  pd = array(dim= c(gCount,4) )
  for (j in 1:gCount){
    # '... Calculate probability density on the Z scale
    tmp <- KernelPD(Zi, -1, Zgrain, Zerr, Num)
    Dsum <- tmp[1]
    SEDsum <- tmp[2]
    Dsum = 100 * barWidth * Dsum / Num
    SEDsum = 100 * barWidth * SEDsum / Num
    pd[j, 1] = TaufromZ(Zi)
    pd[j, 2] = Dsum
    pd[j, 3] = Dsum - SEDsum
    pd[j, 4] = Dsum + SEDsum
    atZ[j]  = Zi
    
    Zi = Zi + zWidth
  }
  
  dat2 = data.frame(atZ, pd[,1],pd[,2],pd[,3],pd[,4],log(pd[,1]))
  colnames(dat2)<-c("Zi","age", "mean","upper", "lower","logAges")
  
  return( dat2 )
  
}


extractPeakResults <- function ( resultsOutput , deltaBIC){
  peakResults<- c(deltaBIC, resultsOutput$logLike)
  for(i in 1:resultsOutput$PkNum){
    peakResults <- c(peakResults, resultsOutput$PeakAgeResults[i], resultsOutput$PeakAgeCI95min[i] , resultsOutput$PeakAgeCI95plus[i] , resultsOutput$PkFrac[i])
  }
  return(peakResults)
}


getStErrRhoDFromRhoDandND <- function(rhoD, nD){
  stErrRhoD = rhoD/sqrt(nD)
}


# Section 3.8 -  the Ni and Ns in the book are the expected values !!
Chi2ageHeterogeneityTest <- function(FTdataset, benchmarkData=FALSE){
  
  if(benchmarkData){
#     This benchmark data comes from Table 3.2 in the Galbraith book
    nS=c(0,2,18,2,10,3,4,20,52,2,1,6,256,52,3,10,2,7,1,14,15,14,8,22,16,34,14,6,13,127)
    nI=c(11,11,28,4,78,22,8,57,129,7,9,16,220,134,11,17,5,23,10,43,44,25,28,69,29,51,56,9,22,213)
#    Quoted solution: Chi2=152.0 ;  df = 29  ;  p-Value = <<< 0.01  => single trueAge
  } else {
    nI = FTdataset$nI
    nS = FTdataset$nS
  }
  
  df = length(nI)-1

  nImean = mean(nI)
  nSmean = mean(nS)
  
  Chi2=0
  for(i in 1:length(nI)){
    Chi2 = Chi2 + (nS[i]*nImean - nI[i]*nSmean)**2 / (nI[i]+nS[i]) 
  }
  Chi2 = Chi2/(nImean*nSmean)
  pValue = pchisq(Chi2, df, lower.tail=FALSE)
  
  if(pValue<0.01){
    comment="pValue<0.01 : Strong evidence against a common age model"
  } else if (pValue<0.05) {
    comment="pValue<0.05 : Moderate evidence against a common age model"
  } else {
    comment = "pValue>0.05 : Data consistent with a common age model"
  }

  tmp = list(Chi2=Chi2 , df=df , pValue=pValue, significant=comment)
  return(tmp)
}


