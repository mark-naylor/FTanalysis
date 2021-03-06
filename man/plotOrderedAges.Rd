\name{plotOrderedAges}
\alias{plotOrderedAges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
New plotting function to plot ordered grain ages with uncertainty
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotOrderedAges(FTdataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
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
function (FTdataset) 
{
    nS = FTdataset$nS
    nI = FTdataset$nI
    zeta = FTdataset$Zeta
    rhoD = FTdataset$rhoD
    nGrains = FTdataset$nGrain
    grainAge = ((1/1.55125e-10) * log((1.55125e-10 * (nS/nI) * 
        rhoD * zeta * 0.5) + 1))/1e+06
    n = order(grainAge)
    index = seq(1, nGrains)
    df = data.frame(ages = grainAge[n], index = index, lowerAge = FTdataset$grainAges95min[n], 
        upperAge = FTdataset$grainAges95max[n])
    ggplot(df) + geom_point(aes(x = ages, y = index)) + geom_errorbarh(aes(x = ages, 
        xmax = upperAge, xmin = lowerAge, y = index)) + xlim(0, 
        100)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
