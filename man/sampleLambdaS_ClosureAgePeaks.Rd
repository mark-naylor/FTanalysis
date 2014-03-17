\name{sampleLambdaS_ClosureAgePeaks}
\alias{sampleLambdaS_ClosureAgePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleLambdaS_ClosureAgePeaks(nI, closureAges, b, probs = NULL, nAges = NULL, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nI}{
%%     ~~Describe \code{nI} here~~
}
  \item{closureAges}{
%%     ~~Describe \code{closureAges} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{probs}{
%%     ~~Describe \code{probs} here~~
}
  \item{nAges}{
%%     ~~Describe \code{nAges} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (nI, closureAges, b, probs = NULL, nAges = NULL, plot = FALSE) 
{
    nGrains = length(nI)
    if (is.null(nAges)) {
        nAges = length(closureAges)
    }
    if (is.null(probs)) {
        probs = rep(1/nAges, nAges)
    }
    whichDist = sample(1:nAges, nGrains, replace = T, prob = probs)
    lambda_nS_true = nI * (exp(lambda.U238 * closureAges[whichDist] * 
        1e+06) - 1)/b
    modelType = paste(nAges, " peak synthetic", sep = "")
    print(modelType)
    if (plot) {
        df3 = data.frame(nI = nI, nS = lambda_nS_true, whichDist = whichDist)
        df3$y1 <- df3$nI * (exp(lambda.U238 * closureAges[1] * 
            1e+06) - 1)/b
        df3$y2 <- df3$nI * (exp(lambda.U238 * closureAges[2] * 
            1e+06) - 1)/b
        synthPlot2 <- ggplot(df3) + geom_point(aes(x = nI, y = nS, 
            color = as.factor(whichDist))) + scale_colour_manual(values = c(`1` = "red", 
            `2` = "blue", `3` = "green")) + geom_line(aes(x = nI, 
            y = y1), color = "red") + geom_line(aes(x = nI, y = y2), 
            color = "blue")
        multiplot(synthPlot2, cols = 2)
    }
    thisSample = list(lambdaS = lambda_nS_true, closureAges = closureAges[whichDist], 
        whichDist = whichDist, modelType = modelType, peaksAt = closureAges)
    return(thisSample)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
