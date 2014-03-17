\name{syntheticSummaryPlot}
\alias{syntheticSummaryPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
syntheticSummaryPlot(FTdataset, synthetic_LambdaS, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{synthetic_LambdaS}{
%%     ~~Describe \code{synthetic_LambdaS} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (FTdataset, synthetic_LambdaS, b) 
{
    if (synthetic_LambdaS$modelType == "Step function synthetic") {
        print("Using step plot")
        df3 = data.frame(nI = FTdataset$nI, nS = FTdataset$nS, 
            closureAges = synthetic_LambdaS$closureAges)
        df3$y1 <- df3$nI * (exp(lambda.U238 * synthetic_LambdaS$minAge * 
            1e+06) - 1)/b
        df3$y2 <- df3$nI * (exp(lambda.U238 * synthetic_LambdaS$maxAge * 
            1e+06) - 1)/b
        syntheticSummaryPlot <- ggplot(df3) + geom_point(aes(x = nI, 
            y = nS, color = closureAges)) + geom_line(aes(x = nI, 
            y = y1), color = "red") + geom_line(aes(x = nI, y = y2), 
            color = "blue")
    }
    else {
        print("Using peak plot")
        df3 = data.frame(nI = FTdataset$nI, nS = FTdataset$nS, 
            whichDist = synthetic_LambdaS$whichDist, trueAge = synthetic_LambdaS$closureAges)
        df3$y1 <- df3$nI * (exp(lambda.U238 * synthetic_LambdaS$peaksAt[1] * 
            1e+06) - 1)/b
        df3$y2 <- df3$nI * (exp(lambda.U238 * synthetic_LambdaS$peaksAt[2] * 
            1e+06) - 1)/b
        df3$y1_C975 <- qpois(0.975, df3$y1)
        df3$y1_C025 <- qpois(0.025, df3$y1)
        df3$y2_C975 <- qpois(0.975, df3$y2)
        df3$y2_C025 <- qpois(0.025, df3$y2)
        df3$y1Color <- df3$nI/df3$nI * min(synthetic_LambdaS$closureAges)
        df3$y2Color <- df3$nI/df3$nI * max(synthetic_LambdaS$closureAges)
        syntheticSummaryPlot <- ggplot(df3) + geom_point(aes(x = nI, 
            y = nS, color = as.factor(trueAge))) + geom_path(aes(x = nI, 
            y = y1, colour = as.factor(y1Color)), linetype = "dashed") + 
            geom_path(aes(x = nI, y = y2, colour = as.factor(y2Color)), 
                linetype = "dashed") + ggtitle("Plot of synthetic dataset") + 
            scale_colour_discrete(guide = guide_legend(title = "Closure Age")) + 
            theme(legend.position = "bottom")
    }
    return(syntheticSummaryPlot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
