makeFTdataset <- function(nS, nI, nD=NULL, rhoD, relErrRhoD, c, K, Zeta, relErrZeta, SqSize = NULL, 	geomFactor = 0.5){
  
  if(length(nI)==length(nS)) { 
    nGrain = length(nS) 
  } else { 
    print( "ERROR: There are different numbers of nI and nS")
  }
  
  NNi = K * nI
  nT = nS + NNi
  
  # '... Calculate Z and grain-only standard error
  NNs = nS   
  NNi = NNi 
  b <- Zeta*rhoD*0.5*LamdaD
  Zgrain = log(b * NNs / NNi)
  Zerr = sqrt(1 / NNs + 1 / NNi)
  QMeanZerr = sqrt ( sum( Zerr * Zerr ) / nGrain )
  
  stErrRhoD <- relErrRhoD * rhoD
  nD <- (rhoD/stErrRhoD)**2
  
  grainAges = TaufromZ(Zgrain)
  grainAges95min = TaufromZ(Zgrain+2* Zerr)
  grainAges95max = TaufromZ(Zgrain-2* Zerr)
  
  FTdataset = list(nS=nS, nI=nI, nD=nD, nT=nT, grainAges = grainAges, rhoD=rhoD, relErrRhoD = relErrRhoD, Zeta = Zeta, relErrZeta = relErrZeta, c=c, nGrain=nGrain, Zgrain = Zgrain, Zerr = Zerr, QMeanZerr= QMeanZerr, b=b, nD=nD, grainAges95min= grainAges95min, grainAges95max= grainAges95max)
  
  class(FTdataset) <- "FTdataset"
  return( FTdataset )
}

# plot.FTdataset <-function(FTdataset) {	
# FTdatasetPlot <- 
# return(FTdatasetPlot)
# }

# # # # # # # # # # # # # # # # # # # # # # # # # # #

dataSummaryPlot <- function(FTdataset, b){
  
  df3 = data.frame(nI=FTdataset$nI, nS=FTdataset$nS)
  
  dataSummaryPlot <- ggplot(df3) + 
    geom_point( aes(x=nI, y=nS ) ) +
    ggtitle("Plot of real dataset")
  
  # geom_path( aes(x=nI, y=y1, colour=as.factor(y1Color)), linetype="dashed" )+
  # scale_colour_discrete(guide = guide_legend(title = "Closure Age"))+ 
  # theme(legend.position="bottom")
  
  return(dataSummaryPlot)
  
}

