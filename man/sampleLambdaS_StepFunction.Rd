\name{sampleLambdaS_StepFunction}
\alias{sampleLambdaS_StepFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleLambdaS_StepFunction(nI, minAge, maxAge, b, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nI}{
%%     ~~Describe \code{nI} here~~
}
  \item{minAge}{
%%     ~~Describe \code{minAge} here~~
}
  \item{maxAge}{
%%     ~~Describe \code{maxAge} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (nI, minAge, maxAge, b, plot = FALSE) 
{
    nGrains = length(nI)
    closureAges = runif(nGrains, minAge, maxAge)
    lambda_nS_true = nI * (exp(lambda.U238 * closureAges * 1e+06) - 
        1)/b
    modelType = paste("Step function synthetic", sep = "")
    print(modelType)
    if (plot) {
        df3 = data.frame(nI = nI, nS = lambda_nS_true, closureAges = closureAges)
        df3$y1 <- df3$nI * (exp(lambda.U238 * minAge * 1e+06) - 
            1)/b
        df3$y2 <- df3$nI * (exp(lambda.U238 * maxAge * 1e+06) - 
            1)/b
        synthPlot2 <- ggplot(df3) + geom_point(aes(x = nI, y = nS, 
            color = closureAges)) + geom_line(aes(x = nI, y = y1), 
            color = "red") + geom_line(aes(x = nI, y = y2), color = "blue")
        multiplot(synthPlot2, cols = 2)
    }
    sample = list(lambdaS = lambda_nS_true, closureAges = closureAges, 
        modelType = modelType, minAge = minAge, maxAge = maxAge)
    return(sample)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
