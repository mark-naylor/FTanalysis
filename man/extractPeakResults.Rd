\name{extractPeakResults}
\alias{extractPeakResults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to format peak fitting results}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractPeakResults(resultsOutput, deltaBIC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsOutput}{
%%     ~~Describe \code{resultsOutput} here~~
}
  \item{deltaBIC}{
%%     ~~Describe \code{deltaBIC} here~~
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
function (resultsOutput, deltaBIC) 
{
    peakResults <- c(deltaBIC, resultsOutput$logLike)
    for (i in 1:resultsOutput$PkNum) {
        peakResults <- c(peakResults, resultsOutput$PeakAgeResults[i], 
            resultsOutput$PeakAgeCI95min[i], resultsOutput$PeakAgeCI95plus[i], 
            resultsOutput$PkFrac[i])
    }
    return(peakResults)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
