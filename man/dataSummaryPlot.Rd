\name{dataSummaryPlot}
\alias{dataSummaryPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot and return a graphs of nI vs nS}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dataSummaryPlot(FTdataset, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
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
function (FTdataset, b) 
{
    df3 = data.frame(nI = FTdataset$nI, nS = FTdataset$nS)
    dataSummaryPlot <- ggplot(df3) + geom_point(aes(x = nI, y = nS)) + 
        ggtitle("Plot of real dataset")
    return(dataSummaryPlot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
