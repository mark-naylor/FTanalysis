
getStochasticBimodalFit<- function(nI, peaksAt, input, b, probs=NULL){
  
  nGrains = length(nI)
  
  synthetic_LambdaS <-sampleLambdaS_ClosureAgePeaks(nI, peaksAt, b , plot=FALSE, probs=probs)
  
  #   **Generate a sample of nS**
  nS<-sampleNs(synthetic_LambdaS$lambdaS)
  
  #   Create the FTdataset object and the data summary plots
  FTdataset <- makeFTdataset(nS=nS, nI=nI, rhoD=input$RhoD, relErrRhoD=input$RERhoD, c=1, K=1, Zeta=input$Zeta, relErrZeta=input$REZeta, SqSize = NULL, geomFactor = 0.5)
  
  #   trackCountSummaryPlot <- plotTrackCountSummary(FTdataset)
  #   PDplotNoOverlay  <- PDplot(FTdataset, resultsOutput, plotType=1, zeroNsOffset=0.5)
  #   plot2  <- PDplot(FTdataset, resultsOutput, plotType=6, zeroNsOffset=0.5)
  
  #   layout <- matrix(c(1, 2, 3), nrow = 1, byrow = TRUE)
  #   multiplot( trackCountSummaryPlot, PDplotNoOverlay, plot2, layout=layout)
  
  #   Calculate the central age
  centralAge <- BINOMFIT_CentralAge(FTdataset)
  
  #   Perform the binomfit analysis
  #     Run Binomfit analysis for different numbers of peaks and hold these results in a list.
  resultsList <- list()
  tryCatch( { 
    resultsOutput1 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=1, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput1 <- resultsOutput1
  }, error = function(ex) {
    print("Excluding Results 1")
  } )
  tryCatch( { 
    resultsOutput2 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=2, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput2 <- resultsOutput2
  }, error = function(ex) {
    print("Excluding Results 2")
  } )
  tryCatch( { 
    resultsOutput3 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=3, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput3 <- resultsOutput3
  }, error = function(ex) {
    print("Excluding Results 3")
  } )
  tryCatch( { 
    resultsOutput4 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=4, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput4 <- resultsOutput4
  }, error = function(ex) {
    print("Excluding Results 4")
  } )
  tryCatch( { 
    resultsOutput5 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=5, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput5 <- resultsOutput5
  }, error = function(ex) {
    print("Excluding Results 5")
  } )
  tryCatch( { 
    resultsOutput6 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=6, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput6 <- resultsOutput6
  }, error = function(ex) {
    print("Excluding Results 6")
  } )
  
  #   Model comparison    
  ### BIC model comparison
  
  nResults=length(resultsList)  
  
  logLikeArray = c()
  nAgeArray =c()
  for(i in 1:nResults){ 
    logLikeArray <- c(logLikeArray, resultsList[[i]]$logLike) 
    nAgeArray <- c(nAgeArray, resultsList[[i]]$PkNum )
  }
  
  BICarray <- -2*logLikeArray + (2*nAgeArray-1) * log(FTdataset$nGrain)
  deltaBIC <- BICarray-min(BICarray)
  
  BICdf <- data.frame(nAges=nAgeArray, logLike=logLikeArray, BIC=BICarray, deltaBIC=deltaBIC)
  
  #   xt<-xtable(BICdf)
  #   BICmodelComparisonPlot <- ggplot(BICdf, aes(x = nAges, y = deltaBIC))+ geom_point( color="blue") + geom_line()
  
  #   print(xt, type="html")
  #   BICmodelComparisonPlot
  
  n = which(BICarray==min(BICarray)) 
  favouredBIC_nPeaks <- resultsList[[n]]$PkNum
  favouredBIC_ages   <- resultsList[[n]]$PeakAgeResults
  favouredBIC_proportions   <- resultsList[[n]]$PkFrac
  favouredBIC_C95Plus <- resultsList[[n]]$PeakAgeCI95plus
  favouredBIC_C95Minus <- resultsList[[n]]$PeakAgeCI95min
  
  ### Chi2 model comparison (needs a little tidying up)
  ModelCompProb<-c()  ;  F<-c()
  
  for(i in 1:(nResults-1)){
    tmp=getChi2Comp(resultsList[[i]]$ChiSq, resultsList[[i]]$degFreedom, resultsList[[i+1]]$ChiSq, resultsList[[i+1]]$degFreedom) 
    ModelCompProb[i]<-tmp$P
  }
  
  # df = data.frame(Model1=seq(1, (nResults-1)), Model2=seq(2,nResults), Prob_F_ByChanceAlone=ModelCompProb*100)
  
  #   print( xtable(df), type="html") 
  #   ggplot(df, aes(x = Model2, y = Prob_F_ByChanceAlone))+ geom_point( color="blue") + geom_line()
  
  #   Summary Plots of the binomfit results
  #   ------------------------------------
  #     This plot summarises the data in the first 3 figures, and the models fitted by FTanalysis in the remainder. Pairs of PD plots and Radial Plots are shown for models with upto 6 peaks. The prefered model is indicated using the minimum BIC.  
  
  #   makeBinomfitSummaryPlot_4AgeModels(FTdataset, resultsOutput1, resultsOutput2, resultsOutput3, resultsOutput4, ageLabels, dataTrasformStyle="arcsinTransformation")
  # makeBinomfitSummaryPlot_6AgeModels(FTdataset, resultsOutput1, resultsOutput2, resultsOutput3, resultsOutput4, resultsOutput5, resultsOutput6, ageLabels, dataTrasformStyle="arcsinTransformation")
  
  #   Compare FTanalysis results with the Bimodal Ages (in blue)
  #   **Bimodal Ages:**        The input ages are at `r peaksAt`
  #   **FTanalysis:**          `r favouredBIC_nPeaks` peaks fitted at `r sort(favouredBIC_ages[1:favouredBIC_nPeaks])`
  
  #   comparisonPlot <- compareModelSolutions(resultsList, favouredBIC_nPeaks)
  #   comparisonPlot <- comparisonPlot + geom_vline(xintercept=peaksAt, colour="blue") + xlim(0,60)
  
  results = list(peaksAt=peaksAt, nGrains=nGrains, nPeaks=favouredBIC_nPeaks, ages=favouredBIC_ages, proportions=favouredBIC_proportions, C95Plus=favouredBIC_C95Plus, C95Minus=favouredBIC_C95Minus, centralAge=centralAge, BICmodelChoice=BICdf )
  return(results)  
  
}






getStochasticUniformFit<- function(nI, limits, input, b){
  
  nGrains = length(nI)
  
  synthetic_LambdaS <- sampleLambdaS_StepFunction(nI, minAge=min(limits), maxAge=max(limits), b )
  
  #   **Generate a sample of nS**
  nS<-sampleNs(synthetic_LambdaS$lambdaS)
  
  #   **The values of nS are:** `r nS`
  #   **The values of nI are:** `r nI`
  
  #   Create the FTdataset object and the data summary plots
  FTdataset <- makeFTdataset(nS=nS, nI=nI, rhoD=input$RhoD, relErrRhoD=input$RERhoD, c=1, K=1, Zeta=input$Zeta, relErrZeta=input$REZeta, SqSize = NULL, geomFactor = 0.5)
  
  #   trackCountSummaryPlot <- plotTrackCountSummary(FTdataset)
  #   PDplotNoOverlay  <- PDplot(FTdataset, resultsOutput, plotType=1, zeroNsOffset=0.5)
  #   plot2  <- PDplot(FTdataset, resultsOutput, plotType=6, zeroNsOffset=0.5)  
  #   layout <- matrix(c(1, 2, 3), nrow = 1, byrow = TRUE)
  #   multiplot( trackCountSummaryPlot, PDplotNoOverlay, plot2, layout=layout)
  
  #   Calculate the central age
  centralAge <- BINOMFIT_CentralAge(FTdataset)
  
  #   Perform the binomfit analysis
  #   Run Binomfit analysis for different numbers of peaks and hold these results in a list.
  resultsList <- list()
  tryCatch( { 
    resultsOutput1 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=1, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput1 <- resultsOutput1
  }, error = function(ex) {
    print("Excluding Results 1")
  } )
  tryCatch( { 
    resultsOutput2 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=2, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput2 <- resultsOutput2
  }, error = function(ex) {
    print("Excluding Results 2")
  } )
  tryCatch( { 
    resultsOutput3 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=3, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput3 <- resultsOutput3
  }, error = function(ex) {
    print("Excluding Results 3")
  } )
  tryCatch( { 
    resultsOutput4 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=4, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput4 <- resultsOutput4
  }, error = function(ex) {
    print("Excluding Results 4")
  } )
  tryCatch( { 
    resultsOutput5 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=5, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput5 <- resultsOutput5
  }, error = function(ex) {
    print("Excluding Results 5")
  } )
  tryCatch( { 
    resultsOutput6 <- BINOMFIT(FTdataset, peakAgeModel=2, PkNum=6, K=input$K, details=input$details, verbose=FALSE)
    resultsList$resultsOutput6 <- resultsOutput6
  }, error = function(ex) {
    print("Excluding Results 6")
  } )
  
  #   Model comparison
  
  ### BIC model comparison
  nResults=length(resultsList)  
  
  logLikeArray = c()
  nAgeArray =c()
  for(i in 1:nResults){ 
    logLikeArray <- c(logLikeArray, resultsList[[i]]$logLike) 
    nAgeArray <- c(nAgeArray, resultsList[[i]]$PkNum )
  }
  
  BICarray <- -2*logLikeArray + (2*nAgeArray-1) * log(FTdataset$nGrain)
  deltaBIC <- BICarray-min(BICarray)
  
  BICdf <- data.frame(nAges=nAgeArray, logLike=logLikeArray, BIC=BICarray, deltaBIC=deltaBIC)
  
  n = which(BICarray==min(BICarray)) 
  favouredBIC_nPeaks <- resultsList[[n]]$PkNum
  favouredBIC_ages   <- resultsList[[n]]$PeakAgeResults
  favouredBIC_proportions   <- resultsList[[n]]$PkFrac
  favouredBIC_C95Plus <- resultsList[[n]]$PeakAgeCI95plus
  favouredBIC_C95Minus <- resultsList[[n]]$PeakAgeCI95min
  
  ### Chi2 model comparison (needs a little tidying up)
  ModelCompProb<-c()  ;  F<-c()
  
  for(i in 1:(nResults-1)){
    tmp=getChi2Comp(resultsList[[i]]$ChiSq, resultsList[[i]]$degFreedom, resultsList[[i+1]]$ChiSq, resultsList[[i+1]]$degFreedom) 
    ModelCompProb[i]<-tmp$P
  }
  
  df = data.frame(Model1=seq(1, (nResults-1)), Model2=seq(2,nResults), Prob_F_ByChanceAlone=ModelCompProb*100)
  
  #   print( xtable(df), type="html") 
  #   ggplot(df, aes(x = Model2, y = Prob_F_ByChanceAlone))+ geom_point( color="blue") + geom_line()
  
  #   Summary Plots of the binomfit results
  #   ------------------------------------
  #     This plot summarises the data in the first 3 figures, and the models fitted by FTanalysis in the remainder. Pairs of PD plots and Radial Plots are shown for models with upto 6 peaks. The prefered model is indicated using the minimum BIC.  
  
  #   makeBinomfitSummaryPlot_4AgeModels(FTdataset, resultsOutput1, resultsOutput2, resultsOutput3, resultsOutput4, ageLabels, dataTrasformStyle="arcsinTransformation")
  # makeBinomfitSummaryPlot_6AgeModels(FTdataset, resultsOutput1, resultsOutput2, resultsOutput3, resultsOutput4, resultsOutput5, resultsOutput6, ageLabels, dataTrasformStyle="arcsinTransformation")
  
  #   Compare FTanalysis results with the Bimodal Ages (in blue)
  #   **Bimodal Ages:**        The input ages are at `r peaksAt`
  #   **FTanalysis:**          `r favouredBIC_nPeaks` peaks fitted at `r sort(favouredBIC_ages[1:favouredBIC_nPeaks])`
  
  #   comparisonPlot <- compareModelSolutions(resultsList, favouredBIC_nPeaks)
  #   comparisonPlot <- comparisonPlot + geom_vline(xintercept=peaksAt, colour="blue") + xlim(0,60)
  
  results = list(peaksAt=peaksAt, nGrains=nGrains, nPeaks=favouredBIC_nPeaks, ages=favouredBIC_ages, proportions=favouredBIC_proportions, C95Plus=favouredBIC_C95Plus, C95Minus=favouredBIC_C95Minus, centralAge=centralAge, BICmodelChoice=BICdf  )
  return(results)  
  
}