LamdaD <- 1.55125E-10  # 'Total decay constant for 238U (yr^-1)
colourList<-c("purple","green","yellow","red")
colourList<-c("A","B","C","D")

radialPlot <- function(resultsOutput, FTdataset, colourBy=1, ageLabels=c(1,5,10,20,40,80,160), noColor=FALSE, addLegend=FALSE, style="logTransformation", colorAsFactor=TRUE){
  
  PkNum = resultsOutput$PkNum
  nS <- FTdataset$nS
  nI <- FTdataset$nI
  zeta <- FTdataset$Zeta
  rhoD <-FTdataset$rhoD
  z <- FTdataset$Zgrain
  sigma <- FTdataset$Zerr
  Tau <- resultsOutput$PeakAgeResults[1:PkNum]		
  
  if(length(colourBy)==1){ pointColour <- rep(FTdataset$nGrain, colourBy) 
  } else { pointColour <- colourBy }
  
  if( min(FTdataset$nI, FTdataset$nS) == 0 ){  
    print("Radial Plot WARNING: Zero track counts in dataset - Using arcsin transformation")
    style="arcsinTransformation" 
  }
  
  if(style=="logTransformation"){
    z0=sum(z/(sigma* sigma))/sum(1/(sigma* sigma))
    
    zAges <- ZfromTau(Tau)
    
    lines <- zAges-z0
    xMax<-2*c(6*cos(atan(lines)))
    
    zLabels <- ZfromTau(ageLabels)
    
  } else if(style=="arcsinTransformation"){		 		
    z <- asin( sqrt( (nS+0.375)/(nS+nI+0.75) )  )
    sigma <- 0.5* sqrt( 1/(nS+nI+0.5) )
    z0 <- atan( sqrt( sum(nS)/(sum(nI) ) ) )
    
    zAges <- atan( sqrt(  (exp(LamdaD * Tau * 1000000) -1 ) / (0.5 * LamdaD * zeta * rhoD))) 
    
    lines <- zAges-z0
    xMax<-3*2*c(6*cos(atan(lines)))
    
    zLabels <-atan( sqrt(  (exp(LamdaD * ageLabels * 1000000) -1 ) / (0.5 * LamdaD * zeta * rhoD))) 		
  }	
  
  yMax<-c(xMax*lines)
  xCoords <- cbind( xMin=rep(0,PkNum), xMax)
  yCoords <- cbind( yMin=rep(0,PkNum), yMax)
  
  # x1=c(0,xMax[1])
  # y1=c(0,yMax[1])
  # x2=c(0,xMax[2])
  # y2=c(0,yMax[2])
  # x3=c(0,xMax[3])
  # y3=c(0,yMax[3])
  
  x <- 1/sigma
  y <- (z-z0)/sigma
  radius <- (max(sqrt(x*x + y*y))+0.5 )
  xTic <- radius/ sqrt( 1+(zLabels-z0)**2 )
  xTic2 <- xTic + 0.1
  xTic3 <- xTic + 0.2
  yTic <- xTic*(zLabels-z0)
  yTic2 <- xTic2*(zLabels-z0)
  yTic3 <- xTic3*(zLabels-z0)
  tics<-data.frame(xStart=xTic, xEnd=xTic2, yStart=yTic, yEnd=yTic2, xLab=xTic3, yLab=yTic3, label=ageLabels)
  
  yLimMax <- max(y)+2  ;  		yLimMin <- min(y)-2
  xLimMax <- max(x)*1.15
  
  radialDF <- data.frame(x, y, pointColour)
  
  dat <- circleFun(c(0,0), radius, npoints = 100)
  
  radialPlot <- ggplot(radialDF) 
  
  if(noColor){ 
    radialPlot <- radialPlot + geom_point( aes(x=x, y=y) )
  } else { 
    if(colorAsFactor){ 
      radialPlot <- radialPlot + geom_point( aes(x=x, y=y, color=as.factor(pointColour)) )
    } else { 
      radialPlot <- radialPlot + geom_point( aes(x=x, y=y, color=pointColour) ) 
    }
  }
  
  radialPlot <- radialPlot  + geom_path(data=dat, aes(x,y) ) + geom_hline( aes(y=0), linetype="dotted") +
    geom_segment(data=tics, aes(x=xStart, xend=xEnd, y=yStart, yend=yEnd))+
    geom_text(data=tics, aes(x=xLab, y=yLab, label=label, hjust=-0.2), size=3) +
    ggtitle(paste("Radial Plot: ",PkNum," peaks")) + xlab("Precision") + ylab("Standardised estimate")  + 
    ylim(2*yLimMin, 2*yLimMax) + xlim(0, xLimMax) 
  
  
  for (D in 1:PkNum){
    xCoord = c(0,xMax[D])  ;   yMinCoord = c(-2, yMax[D]-2) ;   yMaxCoord =c(2, yMax[D]+2)
    position= data.frame(xCoord, yMinCoord, yMaxCoord, col=D, slope=lines[D] )
    
    if(D==1){
      setColour="purple"
    } else if (D==2){
      setColour="green"
    } else if (D==3){
      setColour="yellow"
    }else if (D==4){
      setColour="red"
    }else if (D==5){
      setColour="blue"
    }else if (D==6){
      setColour="pink"
    }
    
    radialPlot <- radialPlot + 
      geom_ribbon(data= position, aes(x= xCoord, ymin= yMinCoord, ymax= yMaxCoord) , fill= setColour, colour= setColour, alpha=0.2) +
      geom_abline(data=position, aes(intercept=0, slope= slope ), linetype="dashed", colour= setColour)
  }
  
  if(addLegend){ radialPlot<-radialPlot+
                   scale_colour_discrete(guide = guide_legend(title = "Closure Age"))+ 
                   theme(legend.position="top")
  } else {radialPlot<-radialPlot + theme(legend.position = "none") }
  
  return(radialPlot)
}


circleFun <- function(center = c(0,0), radius = 6, npoints = 100){
  tt <- seq(pi/2,3*pi/2,length.out = npoints)
  xx <- center[1] - radius * cos(tt)
  yy <- center[2] + radius * sin(tt)
  return(data.frame(x = xx, y = yy))
}


# # # # # # # # # # # # #
# # # # # # # # # # # # #
# # # # # # # # # # # # #

PDplot <- function(FTdataset, resultsOutput=NULL, plotType=1, zeroNsOffset=0, ageLabels=c(1,2,3,5,7,10,20,50,100,200)){
  
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
  
  Zmin = min(Zgrain)   ;     Zmax = max(Zgrain)
  Ages = TaufromZ(Zgrain)
  
  Num = length(Zgrain)
  QMeanZerr = sum( Zerr * Zerr )
  QMeanZerr = sqrt(QMeanZerr / (Num))
  
  barWidth=0.1			# 'width of histogram bar in Z units
  zWidth = barWidth / 5	# 'interval width for PD plot in Z units
  
  LwLmt = barWidth * (floor(Zmin / barWidth) - 1)
  UpLmt = barWidth * (1 + floor(Zmax / barWidth))

  dat = data.frame(Ages, Zgrain, log(Ages))
  colnames(dat)<-c("Ages", "Zgrain","logAges")


#######################   
  # Calculate probability density distribution based on the data
  # Start and end at 5*QMeanZerr below Zmin and above Zmax, respectively.
  LwLmt = Zmin - (5 * QMeanZerr)
  UpLmt = Zmax + (5 * QMeanZerr)
  Zi = LwLmt
  atZ=c()
  # '  set above: zWidth = BarWidth / 5
  gCount = floor((UpLmt - LwLmt) / zWidth) + 1
  pd = array(dim= c(gCount,5) )
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
    bottom=pd[,1]*0
    
  
  if(plotType==1){
    
    p <- ggplot( dat, aes(Ages)) + geom_histogram(binwidth= barWidth) +
      scale_x_continuous(trans="log", breaks=ageLabels, expand=c(0.001,0))
    
  } else if (plotType==2) {
    
    p  <- ggplot(data=dat2) + geom_line(aes(x=age, y=mean)) + 
      geom_line(aes(x=age, y=lower), colour="red", linetype="dashed")+ 
      geom_line(aes(x=age, y=upper), linetype="dashed", colour="red") + xlim(0, 50)
    
  } else if (plotType==3){
    
    p <- ggplot( dat, aes(Ages)) + geom_histogram(binwidth=1.) + 
      geom_line(aes(x=age, y=firstPeak), data=dat3, colour="red", linetype="dotdash") + 
      geom_line(aes(x=age, y=secondPeak), data=dat3, colour="purple", linetype="dotdash")+ 
      geom_line(aes(x=age, y=thirdPeak), data=dat3, colour="yellow", linetype="dotdash")+ 
      geom_line(aes(x=age, y=total), data=dat3, colour="red") + xlim(0, 50)
    
  } else if (plotType==4){
    
    p <- ggplot(data=dat3) + geom_line(aes(x=age, y=mean), colour="red", linetype="dotdash") + geom_line(aes(x=age, y=lower), colour="purple", linetype="dotdash")+ geom_line(aes(x=age, y=upper), colour="blue", linetype="dotdash")+ geom_line(aes(x=age, y=total), colour="black") + xlim(0, 50)
    
  } else if (plotType==5){
    
    p <- ggplot( dat, aes(Zgrain)) + geom_histogram(binwidth= barWidth) + geom_line( aes(x=Zi, y=mean), data=dat2 , colour="red")+ geom_line(aes(x=Zi, y=lower), data=dat2 ,colour="red", linetype="dashed") + geom_line(aes(x=Zi, y=upper), data=dat2 ,linetype="dashed", colour="red")
    
  } else if (plotType==6){
    
    lims = ageLabels
    dfLims = data.frame(x=lims)
    
    p <- ggplot( dat, aes(Ages)) + geom_histogram(binwidth= barWidth) + 
      scale_x_continuous(trans="log", breaks= ageLabels, expand=c(0.001,0)) + 
      geom_line( aes(x= age, y=mean), data=dat2 , colour="red") + 
      geom_line(aes(x= age, y=lower), data=dat2 ,colour="red", linetype="dashed") + 
      geom_line(aes(x= age, y=upper), data=dat2 ,linetype="dashed", colour="red") +
      xlab("GrainAge [Ma]") + ylab("Count") + ggtitle("PD Plot: Data density")
    
  } else if (plotType==7){
    #######################    
    #  Create density curves for the peak age model  
    PeakAgeResults = resultsOutput$PeakAgeResults
    PeakAgeProposed = resultsOutput$PeakAgeProposed
    PkSDz = resultsOutput$PkSDz
    PkFrac = resultsOutput$PkFrac
    PkNum = resultsOutput$PkNum
    PkZ = ZfromTau(PeakAgeResults)
    
    pdAges = array(dim= c(gCount,2 + PkNum) )
    pdAges[,1]=pd[,1]
    
    Zi = LwLmt
    for (j in 1:gCount){
      SumPkD = 0
      for (i in 1:PkNum) {
        dev = Zi - PkZ[i]
        S = PkSDz[i]
        if (abs(dev) < 5 * S) {
          PkD = (100 * barWidth * PkFrac[i]) * (C0 / S) * exp(-(dev * dev) / (2 * S * S))
        } else {
          PkD = 0
        }
        SumPkD = SumPkD + PkD
        pdAges[j, 1 + i] = PkD
      }
      pdAges[j, i+2] = SumPkD
      Zi = Zi + zWidth
    }
    
    if(PkNum==1){
      dat3 = data.frame(atZ, pdAges[,1],pdAges[,2],pdAges[,3],log(pdAges[,1]), bottom)
      colnames(dat3)<-c("Zi","age", "firstPeak","total","logAges","bottom")  		
    } else if (PkNum==2) {
      dat3 = data.frame(atZ, pdAges[,1],pdAges[,2],pdAges[,3],pdAges[,4],log(pdAges[,1]), bottom)
      colnames(dat3)<-c("Zi","age", "firstPeak","secondPeak","total","logAges","bottom")		
    } else if (PkNum==3) {
      dat3 = data.frame(atZ, pdAges[,1],pdAges[,2],pdAges[,3],pdAges[,4],pdAges[,5],log(pdAges[,1]), bottom)
      colnames(dat3)<-c("Zi","age", "firstPeak","secondPeak", "thirdPeak","total","logAges","bottom")	
    } else if (PkNum==4) {
      dat3 = data.frame(atZ, pdAges[,1],pdAges[,2],pdAges[,3],pdAges[,4],pdAges[,5],pdAges[,6],log(pdAges[,1]), bottom)
      colnames(dat3)<-c("Zi","age", "firstPeak","secondPeak", "thirdPeak","fourthPeak","total","logAges","bottom")		
    } else if (PkNum==5) {
      dat3 = data.frame(atZ, pdAges[,1],pdAges[,2],pdAges[,3],pdAges[,4],pdAges[,5],pdAges[,6],pdAges[,7],log(pdAges[,1]), bottom)
      colnames(dat3)<-c("Zi","age", "firstPeak","secondPeak", "thirdPeak","fourthPeak","fifthPeak","total","logAges","bottom")  	
    } else if (PkNum==6) {
      dat3 = data.frame(atZ, pdAges[,1],pdAges[,2],pdAges[,3],pdAges[,4],pdAges[,5],pdAges[,6],pdAges[,7],pdAges[,8],log(pdAges[,1]), bottom)
      colnames(dat3)<-c("Zi","age", "firstPeak","secondPeak", "thirdPeak","fourthPeak","fifthPeak","sixthPeak","total","logAges","bottom")  	
    }
    
    lims = c(1,2,3,5,7,10,20,30,50,70,100,120,130,150,170,200)
    dfLims = data.frame(x=lims)
    
    p <- ggplot( dat, aes(Ages)) + geom_histogram(binwidth= barWidth) + 
      scale_x_continuous(trans="log", breaks=ageLabels, expand=c(0.001,0)) + 
      geom_ribbon(aes(x=age, ymax=firstPeak, ymin=bottom), data=dat3, fill="purple", alpha=0.3, colour="purple", linetype="dotdash", size=0.8) + geom_line(aes(x=age, y=firstPeak), data=dat3, colour="purple", linetype="dotdash", size=0.8)
    
    if (PkNum>1){
      p <- p + geom_ribbon(aes(x=age, ymax= secondPeak, ymin=bottom), data=dat3, fill="green", alpha=0.3, colour="green", linetype="dotdash", size=0.8) + geom_line(aes(x=age, y=secondPeak), data=dat3, colour="green", linetype="dotdash", size=0.8)
    }
    
    if (PkNum>2){
      p <- p + geom_ribbon(aes(x=age, ymax= thirdPeak, ymin=bottom), data=dat3, fill="yellow", alpha=0.3, colour="yellow", linetype="dotdash", size=0.8) + geom_line(aes(x=age, y=thirdPeak), data=dat3, colour="yellow", linetype="dotdash", size=0.8)
    }
    
    if (PkNum>3){
      p <- p + geom_ribbon(aes(x=age, ymax= fourthPeak, ymin=bottom), data=dat3, fill="red", alpha=0.3, colour="red", linetype="dotdash", size=0.8) + geom_line(aes(x=age, y= fourthPeak), data=dat3, colour="red", linetype="dotdash", size=0.8)
    }
    if (PkNum>4){
      p <- p + geom_ribbon(aes(x=age, ymax= fifthPeak, ymin=bottom), data=dat3, fill="blue", alpha=0.3, colour="blue", linetype="dotdash", size=0.8) + geom_line(aes(x=age, y= fifthPeak), data=dat3, colour="blue", linetype="dotdash", size=0.8)
    }
    if (PkNum>5){
      p <- p + geom_ribbon(aes(x=age, ymax= sixthPeak, ymin=bottom), data=dat3, fill="pink", alpha=0.3, colour="pink", linetype="dotdash", size=0.8) + geom_line(aes(x=age, y= sixthPeak), data=dat3, colour="pink", linetype="dotdash", size=0.8)
    }    
    
    p <- p+ geom_line(aes(x=age, y=total), data=dat3, colour="red", size=0.9) + 
      xlab("GrainAge [Ma]") + ylab("Count") + ggtitle(paste("PD Plot: ",PkNum," peaks"))
  }
  
  return( p ) 
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


syntheticSummaryPlot <- function(FTdataset, synthetic_LambdaS, b){
  
  if(synthetic_LambdaS$modelType=="Step function synthetic") {
    print("Using step plot")	
    df3 = data.frame(nI=FTdataset$nI, nS=FTdataset$nS, closureAges= synthetic_LambdaS$closureAges)
    df3$y1 <- df3$nI * ( exp( lambda.U238* synthetic_LambdaS$minAge*1.E6 ) - 1 ) / b
    df3$y2 <- df3$nI * ( exp( lambda.U238* synthetic_LambdaS$maxAge*1.E6 ) - 1 ) / b
    
    syntheticSummaryPlot <- ggplot(df3) + 
      geom_point( aes(x=nI, y=nS, color=closureAges ) ) +
      geom_line( aes(x=nI, y=y1), color="red" ) +
      geom_line( aes(x=nI, y=y2), color="blue" )
    
  } else {
    print("Using peak plot")
    
    df3 = data.frame(nI=FTdataset$nI, nS=FTdataset$nS, whichDist= synthetic_LambdaS$whichDist, trueAge= synthetic_LambdaS$closureAges)
    df3$y1 <- df3$nI * ( exp( lambda.U238* synthetic_LambdaS$peaksAt[1]*1.E6 ) - 1 ) / b
    df3$y2 <- df3$nI * ( exp( lambda.U238* synthetic_LambdaS$peaksAt[2]*1.E6 ) - 1 ) / b
    df3$y1_C975 <- qpois(0.975, df3$y1)
    df3$y1_C025 <- qpois(0.025, df3$y1)
    df3$y2_C975 <- qpois(0.975, df3$y2)
    df3$y2_C025 <- qpois(0.025, df3$y2)
    df3$y1Color <- df3$nI/df3$nI*min(synthetic_LambdaS$closureAges)
    df3$y2Color <- df3$nI/df3$nI*max(synthetic_LambdaS$closureAges)
    
    syntheticSummaryPlot <- ggplot(df3) + 
      geom_point( aes(x=nI, y=nS, color=as.factor(trueAge) ) ) +
      geom_path( aes(x=nI, y=y1, colour=as.factor(y1Color)), linetype="dashed" )+
      geom_path( aes(x=nI, y=y2, colour=as.factor(y2Color)), linetype="dashed" )+ 
      ggtitle("Plot of synthetic dataset")+ 
      scale_colour_discrete(guide = guide_legend(title = "Closure Age"))+ 
      theme(legend.position="bottom")
  }
  
  return(syntheticSummaryPlot)
  
}


# # # # # # # # # # # # # # # # # # # # # # # # # #

plot4x2AgeHists <- function(ageList1, ageList2, ageList3, ageList5, preferedAgeModels1, preferedAgeModels2, preferedAgeModels3, preferedAgeModels4){
  
  blankPlot = ggplot(peakDF1) + geom_blank()+ xlim(0, 10) + ylim(0, 100)
  
  peakDF1 = data.frame(age=ageList1)
  ageHist1 <- ggplot(peakDF1)+ geom_bar( binwidth=1 , aes(x=age, fill=as.factor(1))) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax) # defaults to stacking
  
  peakDF2 = data.frame(age=ageList2[,1] , order=ageList2[,2])
  ageHist2 <- ggplot(peakDF2) + geom_bar( binwidth=1, aes(x=age, fill=as.factor(order)) ) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax) # defaults to stacking
  
  peakDF3 = data.frame(age=ageList3[,1] , order=ageList3[,2])
  ageHist3 <- ggplot(peakDF3) + geom_bar( binwidth=1, aes(x=age, fill=as.factor(order)) ) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax)  # defaults to stacking
  
  peakDF4 = data.frame(age=ageList4[,1] , order=ageList4[,2])
  ageHist4 <- ggplot(peakDF4) + geom_bar( binwidth=1, aes(x=age, fill=as.factor(order)) ) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax) # defaults to stacking
  
  if(length(preferedAgeModels1>0)){
    peakDF1 = data.frame(age= preferedAgeModels1)
    chosenAgeHist1 <- ggplot(peakDF1)+ geom_bar( binwidth=1 , aes(x=age, fill=as.factor(1))) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax) # defaults to stacking
  } else {chosenAgeHist1= blankPlot}
  
  
  if(length(preferedAgeModels2>0)){
    peakDF2 = data.frame(age= preferedAgeModels2[,1] , order= preferedAgeModels2[,2])
    chosenAgeHist2 <- ggplot(peakDF2) + geom_bar( binwidth=1, aes(x=age, fill=as.factor(order)) ) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax) # defaults to stacking
  } else {chosenAgeHist2= blankPlot}
  
  if(length(preferedAgeModels3>0)){
    peakDF3 = data.frame(age= preferedAgeModels3[,1] , order= preferedAgeModels3[,2])
    chosenAgeHist3 <- ggplot(peakDF3) + geom_bar( binwidth=1, aes(x=age, fill=as.factor(order)) ) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax)  # defaults to stacking
  } else {chosenAgeHist3= blankPlot}
  
  if(length(preferedAgeModels4>0)){
    peakDF4 = data.frame(age= preferedAgeModels4[,1] , order= preferedAgeModels4[,2])
    chosenAgeHist4 <- ggplot(peakDF4) + geom_bar( binwidth=1, aes(x=age, fill=as.factor(order)) ) + geom_vline(xintercept = xmin) + geom_vline(xintercept = xmax) # defaults to stacking
  } else {chosenAgeHist4= blankPlot}
  
  theme_set(theme_bw(10))
  layout <- matrix(c(1, 2, 3, 4,5,6,7,8), nrow = 2, byrow = TRUE)
  multiplot(ageHist1, ageHist2, ageHist3, ageHist4, chosenAgeHist1, chosenAgeHist2, chosenAgeHist3, chosenAgeHist4 ,layout=layout)
  
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # 

plotOrderedAges <- function(FTdataset){
  nS = FTdataset$nS  ;      nI = FTdataset$nI
  zeta =  FTdataset$Zeta
  rhoD = FTdataset$rhoD
  nGrains = FTdataset$nGrain
  grainAge =((1/0.000000000155125)*log((0.000000000155125*(nS/nI)*rhoD*zeta*0.5)+1))/1000000
  n=order(grainAge)
  index=seq(1, nGrains)
  df = data.frame(ages=grainAge[n], index=index, lowerAge=FTdataset$grainAges95min[n], upperAge=FTdataset$grainAges95max[n])
  ggplot(df)+geom_point( aes(x=ages, y=index)) + geom_errorbarh( aes(x=ages, xmax= upperAge, xmin= lowerAge,y=index)) + xlim(0,100)
}


plotTrackCountSummary <- function(FTdataset){
  df = data.frame(nS = FTdataset$nS , nI = FTdataset$nI)
  countPlot <- ggplot(df)+geom_point( aes(x=nS, y=nI)) 
  return(countPlot)
}

makeBinomfitSummaryPlot_4AgeModels <-function(FTdataset, resultsOutput1, resultsOutput2, resultsOutput3, resultsOutput4, ageLabels, dataTrasformStyle="arcsinTransformation"){
  
  BICmodelComparisonPlot <- ggplot(BICdf, aes(x = nAges, y = deltaBIC)) + geom_point( color="blue") + geom_line()
  
  trackCountSummaryPlot <- plotTrackCountSummary(FTdataset)
  PDplotNoOverlay  <- PDplot(FTdataset, resultsOutput1, plotType=1, zeroNsOffset=0.5)
  plot2  <- PDplot(FTdataset, resultsOutput1, plotType=6, zeroNsOffset=0.5)
  
  cols <- rep(1.1, FTdataset$nGrain)
  
  radialPlotPlot1 <- radialPlot(resultsOutput1, FTdataset, colourBy=cols, noColor=FALSE, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot1  <- PDplot(FTdataset, resultsOutput1, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot2 <- radialPlot(resultsOutput2, FTdataset, colourBy=1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot2  <- PDplot(FTdataset, resultsOutput2, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot3 <- radialPlot(resultsOutput3, FTdataset, colourBy= 1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot3  <- PDplot(FTdataset, resultsOutput3, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot4 <- radialPlot(resultsOutput4, FTdataset, colourBy= 1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot4  <- PDplot(FTdataset, resultsOutput4, plotType=7, zeroNsOffset=0.5)
  
  theme_set(theme_bw(10))
  layout <- matrix(c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12), nrow = 3, byrow = TRUE)
  summaryPlot <- multiplot( trackCountSummaryPlot, PDplotNoOverlay, plot2, BICmodelComparisonPlot, PDplot1, radialPlotPlot1, PDplot2, radialPlotPlot2, PDplot3, radialPlotPlot3, PDplot4, radialPlotPlot4, layout=layout)
  
  return(summaryPlot)
}


makeBinomfitSummaryPlot_6AgeModels <-function(FTdataset, resultsOutput1, resultsOutput2, resultsOutput3, resultsOutput4, resultsOutput5, resultsOutput6, ageLabels, dataTrasformStyle="arcsinTransformation"){
  
  BICmodelComparisonPlot <- ggplot(BICdf, aes(x = nAges, y = deltaBIC)) + geom_point( color="blue") + geom_line()
  
  trackCountSummaryPlot <- plotTrackCountSummary(FTdataset)
  PDplotNoOverlay  <- PDplot(FTdataset, resultsOutput1, plotType=1, zeroNsOffset=0.5)
  plot2  <- PDplot(FTdataset, resultsOutput1, plotType=6, zeroNsOffset=0.5)
  
  cols <- rep(1.1, FTdataset$nGrain)
  
  radialPlotPlot1 <- radialPlot(resultsOutput1, FTdataset, colourBy=cols, noColor=FALSE, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot1  <- PDplot(FTdataset, resultsOutput1, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot2 <- radialPlot(resultsOutput2, FTdataset, colourBy=1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot2  <- PDplot(FTdataset, resultsOutput2, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot3 <- radialPlot(resultsOutput3, FTdataset, colourBy= 1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot3  <- PDplot(FTdataset, resultsOutput3, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot4 <- radialPlot(resultsOutput4, FTdataset, colourBy= 1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot4  <- PDplot(FTdataset, resultsOutput4, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot5 <- radialPlot(resultsOutput5, FTdataset, colourBy= 1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot5  <- PDplot(FTdataset, resultsOutput5, plotType=7, zeroNsOffset=0.5)
  
  radialPlotPlot6 <- radialPlot(resultsOutput6, FTdataset, colourBy= 1, ageLabels=ageLabels, colorAsFactor=TRUE, style=dataTrasformStyle)
  PDplot6  <- PDplot(FTdataset, resultsOutput6, plotType=7, zeroNsOffset=0.5)
  
  theme_set(theme_bw(10))
  layout <- matrix(c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16), nrow = 4, byrow = TRUE)
  summaryPlot <- multiplot( trackCountSummaryPlot, PDplotNoOverlay, plot2, BICmodelComparisonPlot, PDplot1, radialPlotPlot1, PDplot2, radialPlotPlot2, PDplot3, radialPlotPlot3, PDplot4, radialPlotPlot4, PDplot5, radialPlotPlot5, PDplot6, radialPlotPlot6, layout=layout)
  
  return(summaryPlot)
}



compareModelSolutions<- function(benchmarkSolutions, resultsList, favouredBIC_nPeaks){

benchType=rep("Benchmark",benchmarkSolutions$nPeaks)
dfBenchmark = data.frame(ages=benchmarkSolutions$ages , upperAge=benchmarkSolutions$ages+2*benchmarkSolutions$ageSE, lowerAge=benchmarkSolutions$ages-2*benchmarkSolutions$ageSE, type=benchType , col="Benchmark Solution")

FTtypeType=rep("1 Peak",1)
if(favouredBIC_nPeaks==1){col="Favoured Model"} else {col="Other Solution"}
resultsOutput = resultsList[[1]]
dfFTanalysis1 = data.frame(ages=resultsOutput$PeakAgeResults , upperAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95plus, lowerAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95min, type=FTtypeType, col=col)

FTtypeType=rep("2 Peaks",2)
if(favouredBIC_nPeaks==2){col="Favoured Model"} else {col="Other Solution"}
resultsOutput = resultsList[[2]]
dfFTanalysis2 = data.frame(ages=resultsOutput$PeakAgeResults , upperAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95plus, lowerAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95min, type=FTtypeType, col=col)

FTtypeType=rep("3 Peaks",3)
if(favouredBIC_nPeaks==3){col="Favoured Model"} else {col="Other Solution"}
resultsOutput = resultsList[[3]]
dfFTanalysis3 = data.frame(ages=resultsOutput$PeakAgeResults , upperAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95plus, lowerAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95min, type=FTtypeType, col=col)

FTtypeType=rep("4 Peaks",4)
if(favouredBIC_nPeaks==4){col="Favoured Model"} else {col="Other Solution"}
resultsOutput = resultsList[[4]]
dfFTanalysis4 = data.frame(ages=resultsOutput$PeakAgeResults , upperAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95plus, lowerAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95min, type=FTtypeType, col=col)

FTtypeType=rep("5 Peaks",5)
if(favouredBIC_nPeaks==5){col="Favoured Model"} else {col="Other Solution"}
resultsOutput = resultsList[[5]]
dfFTanalysis5 = data.frame(ages=resultsOutput$PeakAgeResults , upperAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95plus, lowerAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95min, type=FTtypeType, col=col)

FTtypeType=rep("6 Peaks",6)
if(favouredBIC_nPeaks==6){col="Favoured Model"} else {col="Other Solution"}
resultsOutput = resultsList[[6]]
dfFTanalysis6 = data.frame(ages=resultsOutput$PeakAgeResults , upperAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95plus, lowerAge=resultsOutput$PeakAgeResults+resultsOutput$PeakAgeCI95min, type=FTtypeType, col=col)

resultsComparisonPlot <- ggplot(dfBenchmark) + geom_point( aes(x=ages, y=type)) + geom_errorbarh( aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) + 
  geom_point( data=dfFTanalysis1, aes(x=ages, y=type))+ geom_errorbarh( data=dfFTanalysis1, aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) +
  geom_point( data=dfFTanalysis2, aes(x=ages, y=type))+ geom_errorbarh( data=dfFTanalysis2, aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) +
  geom_point( data=dfFTanalysis3, aes(x=ages, y=type))+ geom_errorbarh( data=dfFTanalysis3, aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) +
  geom_point( data=dfFTanalysis4, aes(x=ages, y=type))+ geom_errorbarh( data=dfFTanalysis4, aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) +
  geom_point( data=dfFTanalysis5, aes(x=ages, y=type))+ geom_errorbarh( data=dfFTanalysis5, aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) +
  geom_point( data=dfFTanalysis6, aes(x=ages, y=type))+ geom_errorbarh( data=dfFTanalysis6, aes(x=ages, xmax=upperAge, xmin=lowerAge, y=type, colour=col)) 

return(resultsComparisonPlot)
}


