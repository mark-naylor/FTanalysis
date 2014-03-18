\name{radialPlot}
\alias{radialPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot and return a radial plot using ggplot2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
radialPlot(resultsOutput, FTdataset, colourBy = 1, ageLabels = c(1, 5, 10, 20, 40, 80, 160), noColor = FALSE, addLegend = FALSE, style = "logTransformation", colorAsFactor = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsOutput}{
%%     ~~Describe \code{resultsOutput} here~~
}
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{colourBy}{
%%     ~~Describe \code{colourBy} here~~
}
  \item{ageLabels}{
%%     ~~Describe \code{ageLabels} here~~
}
  \item{noColor}{
%%     ~~Describe \code{noColor} here~~
}
  \item{addLegend}{
%%     ~~Describe \code{addLegend} here~~
}
  \item{style}{
%%     ~~Describe \code{style} here~~
}
  \item{colorAsFactor}{
%%     ~~Describe \code{colorAsFactor} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (resultsOutput, FTdataset, colourBy = 1, ageLabels = c(1, 
    5, 10, 20, 40, 80, 160), noColor = FALSE, addLegend = FALSE, 
    style = "logTransformation", colorAsFactor = TRUE) 
{
    PkNum = resultsOutput$PkNum
    nS <- FTdataset$nS
    nI <- FTdataset$nI
    zeta <- FTdataset$Zeta
    rhoD <- FTdataset$rhoD
    z <- FTdataset$Zgrain
    sigma <- FTdataset$Zerr
    Tau <- resultsOutput$PeakAgeResults[1:PkNum]
    if (length(colourBy) == 1) {
        pointColour <- rep(FTdataset$nGrain, colourBy)
    }
    else {
        pointColour <- colourBy
    }
    summary(resultsOutput)
    if (min(FTdataset$nI, FTdataset$nS) == 0) {
        print("Radial Plot WARNING: Zero track counts in dataset - Using arcsin transformation")
        style = "arcsinTransformation"
    }
    if (style == "logTransformation") {
        z0 = sum(z/(sigma * sigma))/sum(1/(sigma * sigma))
        zAges <- ZfromTau(Tau)
        lines <- zAges - z0
        xMax <- 2 * c(6 * cos(atan(lines)))
        zLabels <- ZfromTau(ageLabels)
    }
    else if (style == "arcsinTransformation") {
        z <- asin(sqrt((nS + 0.375)/(nS + nI + 0.75)))
        sigma <- 0.5 * sqrt(1/(nS + nI + 0.5))
        z0 <- atan(sqrt(sum(nS)/(sum(nI))))
        zAges <- atan(sqrt((exp(LamdaD * Tau * 1e+06) - 1)/(0.5 * 
            LamdaD * zeta * rhoD)))
        lines <- zAges - z0
        xMax <- 3 * 2 * c(6 * cos(atan(lines)))
        zLabels <- atan(sqrt((exp(LamdaD * ageLabels * 1e+06) - 
            1)/(0.5 * LamdaD * zeta * rhoD)))
    }
    yMax <- c(xMax * lines)
    xCoords <- cbind(xMin = rep(0, PkNum), xMax)
    yCoords <- cbind(yMin = rep(0, PkNum), yMax)
    x <- 1/sigma
    y <- (z - z0)/sigma
    radius <- (max(sqrt(x * x + y * y)) + 0.5)
    xTic <- radius/sqrt(1 + (zLabels - z0)^2)
    xTic2 <- xTic + 0.1
    xTic3 <- xTic + 0.2
    yTic <- xTic * (zLabels - z0)
    yTic2 <- xTic2 * (zLabels - z0)
    yTic3 <- xTic3 * (zLabels - z0)
    tics <- data.frame(xStart = xTic, xEnd = xTic2, yStart = yTic, 
        yEnd = yTic2, xLab = xTic3, yLab = yTic3, label = ageLabels)
    yLimMax <- max(y) + 2
    yLimMin <- min(y) - 2
    xLimMax <- max(x) * 1.15
    print(x)
    radialDF <- data.frame(x, y, pointColour)
    dat <- circleFun(c(0, 0), radius, npoints = 100)
    radialPlot <- ggplot(radialDF)
    if (noColor) {
        radialPlot <- radialPlot + geom_point(aes(x = x, y = y))
    }
    else {
        if (colorAsFactor) {
            radialPlot <- radialPlot + geom_point(aes(x = x, 
                y = y, color = as.factor(pointColour)))
        }
        else {
            radialPlot <- radialPlot + geom_point(aes(x = x, 
                y = y, color = pointColour))
        }
    }
    radialPlot <- radialPlot + geom_path(data = dat, aes(x, y)) + 
        geom_hline(aes(y = 0), linetype = "dotted") + geom_segment(data = tics, 
        aes(x = xStart, xend = xEnd, y = yStart, yend = yEnd)) + 
        geom_text(data = tics, aes(x = xLab, y = yLab, label = label, 
            hjust = -0.2), size = 3) + ggtitle("Radial Plot") + 
        xlab("Precision") + ylab("Standardised estimate") + ylim(2 * 
        yLimMin, 2 * yLimMax) + xlim(0, xLimMax)
    for (D in 1:PkNum) {
        xCoord = c(0, xMax[D])
        yMinCoord = c(-2, yMax[D] - 2)
        yMaxCoord = c(2, yMax[D] + 2)
        position = data.frame(xCoord, yMinCoord, yMaxCoord, col = D, 
            slope = lines[D])
        if (D == 1) {
            setColour = "purple"
        }
        else if (D == 2) {
            setColour = "green"
        }
        else if (D == 3) {
            setColour = "yellow"
        }
        else if (D == 4) {
            setColour = "red"
        }
        radialPlot <- radialPlot + geom_ribbon(data = position, 
            aes(x = xCoord, ymin = yMinCoord, ymax = yMaxCoord), 
            fill = setColour, colour = setColour, alpha = 0.2) + 
            geom_abline(data = position, aes(intercept = 0, slope = slope), 
                linetype = "dashed", colour = setColour)
    }
    if (addLegend) {
        radialPlot <- radialPlot + scale_colour_discrete(guide = guide_legend(title = "Closure Age")) + 
            theme(legend.position = "top")
    }
    else {
        radialPlot <- radialPlot + theme(legend.position = "none")
    }
    return(radialPlot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
