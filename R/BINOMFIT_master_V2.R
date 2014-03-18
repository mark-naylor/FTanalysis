C0 = .3989423         # 'Used in Gaussian equation, equals 1/sqrt(2*Pi)
LamdaD = 1.55125E-10  # 'Total decay constant for 238U (yr^-1)


# # # # # # Returns Pooled Age +/1 1sd - copied from spreadsheet calculations
# Used when we can assume that all of the grain ages are the same
BINOMFIT_PooledAge <- function(FTdataset, details=NULL, K=1){
  
  nS_Sum <- sum(FTdataset$nS)
  nI_Sum <- sum(FTdataset$nI)
  
  beta <- log(nS_Sum/nI_Sum)
  varBeta <- 1/nS_Sum + 1/nI_Sum	
  
  age <- 1/LamdaD * log(1+ FTdataset$rhoD *FTdataset$Zeta *LamdaD*0.5*nS_Sum/nI_Sum)*1E-6
  seAge <- age * sqrt( 1/nS_Sum + 1/nS_Sum + 1/FTdataset$nD + (FTdataset$relErrZeta/FTdataset$Zeta)**2)
  
  results <- list(age=age, error=seAge)
  return(results)
}


BINOMFIT_CentralAge  <- function (FTdataset) {
  # Read in data
  nT <- FTdataset$nT  ;  nS <- FTdataset$nS  ;  nI <- FTdataset$nI
  
  # Process data
  y <- nS/nT  ;  z <- log10((nS+0.5)/(nI+0.5))
  
  # Initialise arrars
  sigma<-c() ; eta<-c() ; w <- c()  ;  wA <-c()  ;  wB<-c()
  
  # Set initial guesses (Galbraith and Laslett, 1993)
  sigma[1] <- 0.6 * sd(z)  ;  eta[1] <- sum(nS)/sum(nT)
  
  
  w <- nT /(  (  eta[1]*(1-eta[1])  )  +  (nT-1)*(  eta[1]**2  ) * (  (1-eta[1])**2  ) * sigma[1]**2 )
  wA <- w*w*( y - eta[1] )**2
  wB <- w * y
  
  for(i in 2:20){
    eta[i]     <- sum(wB) / sum(w)
    sigma[i] <- sigma[i-1] * sqrt(  sum(wA) / sum(w) )
    w     <-  nT /(  (eta[i]*(1-eta[i]))  +  (nT-1)*(eta[i]* eta[i])*((1-eta[i])*(1-eta[i]))* sigma[i]* sigma[i])
    wA   <-  w**2  * (y-eta[i])**2
    wB   <-  w * y
  }
  
  nSExp =(sum(nS)*(nS+nI))/(sum(nS)+sum(nI))
  nIExp =(sum(nI)*(nS+nI))/(sum(nS)+sum(nI))
  
  ChiFactor=(((nS-nSExp)**2)/nSExp)+(((nI-nIExp)**2)/nIExp)
  ChiSquared = sum(ChiFactor)
  degFreedom = length(nI)-1
  print(degFreedom)
  
  centralAge <- (1/LamdaD) * (log(( LamdaD*(eta[20]/(1-eta[20]))*FTdataset$rhoD*FTdataset$Zeta*0.5)+1))/1000000
  centralAgeStError <-centralAge* sqrt( (1 / ( eta[20]**2 * (1-eta[20])**2 * sum(w))) +(1/FTdataset$nD) + ((FTdataset$relErrZeta)**2) )
  
  dispersion <- 100*sigma[20]
  
  results = list(centralAge = centralAge, centralAgeStError =centralAgeStError , ChiSquared=ChiSquared, degFreedom=degFreedom, eta=eta, dispersion=dispersion)
  return( results )
  
}



BINOMFIT <- function (FTdataset, TrialAges, PkNum, details=NULL, K=1){
  
  PeakAge = array( dim=( 10 ) )
  SEPeakAge = array( dim=c(10, 4))
  PkTheta = array( dim=(10))
  PkFrac = array( dim=(10))
  SEPkfrac = array( dim=(10))
  PkSDz = array( dim=(10) )
  
  RhoD = FTdataset$rhoD ; relErrRhoD = FTdataset$relErrRhoD
  Zeta = FTdataset$Zeta ; relErrZeta = FTdataset$relErrZeta
  Zgrain = FTdataset$Zgrain  ;  Zerr = FTdataset$Zerr
  Ns = FTdataset$nS  ;  Ni = FTdataset$nI  ;  Nt = FTdataset$nT
  grainAges = FTdataset$grainAges
  QMeanZerr = FTdataset$QMeanZerr
  nGrain = FTdataset$nGrain
  
  peak = array(dim=c(PkNum,3))
  
  
  GF = .5               # 'Geometry factor
  KernelFactor = .6     # 'Scaling factor used in PD plots
  PkWidthFactor = sqrt(1 + KernelFactor * KernelFactor)
  
  Label1 = "=================BinomFit Program  v. 1.8  (Brandon 9/15/04)=================="
  print(  Label1 )
  print(  "This program determines best-fit peaks based on a binomial model described by" )
  print(  "Galbraith and Green (1990, Nuclear Tracks and Radiation Measurements, v. 17," )
  print(  "p. 197-206). The user must give an initial guess for the number of peaks" )
  print(  "and the approximate ages of the peaks. The fit is found using the EM method" )
  print(  "outlined in Galbraith (1988, Technometrics, v. 30, p. 280)." )
  print( " " )
  print(  "Note that BINOMFIT does not include the routine used in ZETAAGE to " )
  print(  "calculate  precise ages for grains with low track densities. As a result," )
  print(  "the range in grain ages reported here may differ somewhat from that reported" )
  print(  "in ZETAAGE but this difference has no influence on the BINOMFIT results." )
  print( " " )
  print(  "Methods and fission track constants are after Hurford and Green (1983)." )
  print(  paste( "Total decay constant for 238U (yr^-1) = ", LamdaD ) )
  print(  paste( "Geometry factor = ", GF ) )
  print( "" )
  
  Zmax = max(Zgrain)
  Zmin = min(Zgrain)
  
  grainAgeMin = min(grainAges)
  grainAgeMax = max(grainAges)
  
  EstPkWidth = PkWidthFactor * QMeanZerr
  
  # '... The PKFind routine gets the initial guess for the peak parameters.
  # ' PkNum% = number of peaks.
  # ' PkTheta(1 to PkNum%) = Theta ratios for peak ages
  # ' Pkfrac(1 to PkNum%) = Proportion of total distribution for each peak
  
  # # # # SEE PAGE 88 for alternative initial conditions
  
  PeakAge[1:PkNum]=TrialAges[1:PkNum]
  
  # outputFileName = paste(dir,filename,"_nAges",PkNum, sep="")
  
  PkTheta = Theta(PeakAge, Zeta, GF, RhoD)
  
  # '... Make a guess of proportional size of each peak
  for (I in 1:PkNum){
    # '... Calculate probability density at the peak.
    Zpk = ZfromTau( PeakAge[I] )
    PkDen = SumPD(Zpk, 0, Zgrain, Zerr, nGrain, KernelFactor)
    # '... Fraction of total grains = peak density * std dev *sqrt(2*Pi)/ nGrain%
    # ' Use a SD(z) = 0.15 to estimate the proportion of each peak
    # ' relative to the total distribution.
    PkFrac[I] = PkDen * .15 / (C0 * nGrain)
    # '... Set a lower limit for PkFrac(I%) to avoid division-by-zero problems.
    if (PkFrac[I] < .001) { PkFrac[I] = .001 }
  }
  
  
  
  print( " " )
  
  # '... Begin routine to find best-fit peaks
  print(   Label1 )
  print(   "FILE HEADER:" )
  print(   details )
  print( "" )
  print(   "FIT OPTION: Best-fit peaks using the binomial model of Galbraith and Green" )
  print(  "" )
  print(   "---------------------INITIAL GUESS for MODEL PARAMETERS------------------------" )
  print(  paste( "NUMBER OF PEAKS TO FIT = ", PkNum ) )
  print(  "PEAK #)  PEAK AGE     THETA    FRACTION(%)  COUNT" )
  
  for (I in 1:PkNum) {print(  paste( I, PeakAge[I], PkTheta[I], PkFrac[I] * 100, nGrain * PkFrac[I])) }
  print(   paste( "                    TOTAL RANGE for GRAIN AGES = ", grainAgeMin, " to ", grainAgeMax, " Ma" ))
  print(   paste( "                              NUMBER OF GRAINS = ", nGrain ))
  print(   paste( "                    DEGREES OF FREEDOM for FIT = ", nGrain - (2 * PkNum - 1) ))
  print(   paste( "         AVERAGE OF THE SE(Z)'S for THE GRAINS = ", QMeanZerr ))
  print(   paste( "ESTIMATED WIDTH OF PEAKS IN PD PLOT IN Z UNITS = ", EstPkWidth ))
  
  
  
  # '... Find best-fit peaks using EM method
  print(  "   *** EM routine is finding a set of best-fit peaks ***" )
  results <- EM(PkTheta, PkFrac, PkNum, nGrain, Ns, Ni, Nt)
  
  
  
  # '... Calculate standard errors
  results2 <- StdErr(FTdataset=FTdataset, PkTheta=results$PkTheta, PkFrac=results$PkFrac, PkNum=PkNum, GF=GF)
  # FTdataset=FTdataset; PkTheta=results$PkTheta; PkFrac=results$PkFrac; PkNum=PkNum; GF=GF	
  
  # '... Output results
  print(  "" )
  print(   "------------------------PARAMETERS for BEST-FIT PEAKS--------------------------" )
  print(  "{Standard error for peak age includes group error}" )
  print(  paste( "{Peak width is for PD plot assuming a kernel factor = ", KernelFactor, "}"))
  print(   "PEAK   ----PEAK AGE & CONF. INTERVAL (MA)---- PEAK WIDTH --GRAINS IN PEAK----")
  print(  "NUMBER   MEAN (----68% CI---) (----95% CI---)    W(Z)    FRAC(%) SE(%)  COUNT")
  TotalNum = 0
  TotalFrac = 0
  D1 <- array(dim=PkNum) ; D2 <- array(dim=PkNum) ; D3 <- array(dim=PkNum)
  for (I in 1: PkNum){
    TotalNum = TotalNum + results$PkFrac[I] * nGrain
    TotalFrac = TotalFrac + results$PkFrac[I]
    D1[I] = 100 * results$PkFrac[I]
    D2[I] = 100 * results2$SEPkfrac[I]
    D3[I] = results$PkFrac[I] * nGrain
    PkSDz[I] = PkWidthFactor * results2$PkSDz[I]
    print(  paste( I, results2$PeakAge[I], results2$PeakAgeCI65min[I], results2$PeakAgeCI65plus[I], results2$PeakAgeCI95min[I], results2$PeakAgeCI95plus[I], PkSDz[I], D1[I], D2[I], D3[I], sep="   "))
  }
  print(  paste("Total:",100 * TotalFrac, TotalNum , sep="   ") )
  print(  " ")
  
  
  ChiSq = Chi2(results$PkTheta, results$PkFrac, PkNum, nGrain, Nt, Ns)
  RChiSq = ChiSq / (nGrain - (2 * PkNum - 1))
  
  print(  paste( "        LOG-LIKELIHOOD for BEST FIT = ", results$LogLike[1] ))
  print(  paste("     CHI-SQUARED VALUE for BEST FIT = ", ChiSq ))
  print(  paste("          REDUCED CHI-SQUARED VALUE = ", RChiSq ))
  print(  paste("         DEGREES OF FREEDOM for FIT = ", nGrain - (2 * PkNum - 1) ))
  # print(  paste("  CONDITION NUMBER for COVAR MATRIX = ", Cnum ))
  print(  paste("               NUMBER OF ITERATIONS = ", results$Iter[1] ))
  
  
  
  resultsOutput <- list( PeakAgeResults=results2$PeakAge, PeakAgeCI65min = results2$PeakAgeCI65min, PeakAgeCI65plus = results2$PeakAgeCI65plus, PeakAgeCI95min = results2$PeakAgeCI95min, PeakAgeCI95plus = results2$PeakAgeCI95plus, PkSDz = PkSDz, D1=D1, D2=D2, D3=D3, PeakAgeProposed = PeakAge, PkSDz = PkSDz, PkFrac = results$PkFrac, PkNum = PkNum, PkZ = results$PkZ, logLike=results$LogLike[1], ChiSq=ChiSq, ReducedChiSq=RChiSq, degFreedom = (nGrain - (2 * PkNum - 1) ) )
  
  return(resultsOutput)
}

