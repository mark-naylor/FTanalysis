
lambda.U238 <- 1.55125e-10

sampleLambdaS_ClosureAgePeaks <- function(nI, closureAges, b, probs = NULL, nAges = NULL, plot=FALSE){
  nGrains = length(nI)
  if( is.null(nAges) ){ nAges = length(closureAges) }
  
  if( is.null(probs) ){ probs = rep(1/nAges, nAges) }
  
  whichDist = sample(1:nAges, nGrains, replace=T, prob=probs)
  lambda_nS_true = nI * ( exp( lambda.U238* closureAges[whichDist]*1.E6 ) - 1 ) / b 
  
  modelType=paste(nAges, " peak synthetic", sep="")
  print(modelType)
  
  if(plot){
    df3 = data.frame(nI=nI, nS=lambda_nS_true, whichDist=whichDist)
    df3$y1 <- df3$nI * ( exp( lambda.U238* closureAges[1]*1.E6 ) - 1 ) / b
    df3$y2 <- df3$nI * ( exp( lambda.U238* closureAges[2]*1.E6 ) - 1 ) / b
    
    synthPlot2 <- ggplot(df3) + 
      geom_point( aes(x=nI, y=nS, color=as.factor(whichDist) ) ) +
      scale_colour_manual( values = c("1" = "red","2" = "blue","3" = "green") ) + 
      geom_line( aes(x=nI, y=y1), color="red" ) +
      geom_line( aes(x=nI, y=y2), color="blue" )
    multiplot( synthPlot2, cols=2)
  }
  
  thisSample = list(lambdaS = lambda_nS_true, closureAges = closureAges[whichDist], whichDist = whichDist, modelType = modelType, peaksAt= closureAges)
  return(thisSample)
}

# Function to generate a synthetic set of average spontaneous track counts from a uniform age distribution
# between a  minimum and maximum age
sampleLambdaS_StepFunction <- function(nI, minAge, maxAge, b, plot=FALSE){
  nGrains = length(nI)
  
  closureAges = runif(nGrains, minAge, maxAge)
  
  lambda_nS_true = nI * ( exp( lambda.U238* closureAges*1.E6 ) - 1 ) / b 
  
  modelType=paste("Step function synthetic", sep="")
  print(modelType)
  
  if(plot){
    df3 = data.frame(nI=nI, nS=lambda_nS_true, closureAges=closureAges)
    df3$y1 <- df3$nI * ( exp( lambda.U238* minAge*1.E6 ) - 1 ) / b
    df3$y2 <- df3$nI * ( exp( lambda.U238* maxAge*1.E6 ) - 1 ) / b
    
    synthPlot2 <- ggplot(df3) + 
      geom_point( aes(x=nI, y=nS, color=closureAges ) ) +
      geom_line( aes(x=nI, y=y1), color="red" ) +
      geom_line( aes(x=nI, y=y2), color="blue" )
    multiplot( synthPlot2, cols=2)
  }
  
  sample = list(lambdaS = lambda_nS_true, closureAges = closureAges, modelType = modelType, minAge=minAge, maxAge=maxAge)
  return(sample)
}


sampleNs <- function( lambdaS ){
  nS <- rpois(length(lambdaS), lambdaS)
  return(nS)
}
