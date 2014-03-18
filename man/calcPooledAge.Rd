\name{calcPooledAge}
\alias{calcPooledAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the pooled age for a set of grains
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcPooledAge(nS, nI, Zeta, ZetaStErr, RhoD, nD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nS}{
%%     ~~Describe \code{nS} here~~
}
  \item{nI}{
%%     ~~Describe \code{nI} here~~
}
  \item{Zeta}{
%%     ~~Describe \code{Zeta} here~~
}
  \item{ZetaStErr}{
%%     ~~Describe \code{ZetaStErr} here~~
}
  \item{RhoD}{
%%     ~~Describe \code{RhoD} here~~
}
  \item{nD}{
%%     ~~Describe \code{nD} here~~
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
function (nS, nI, Zeta, ZetaStErr, RhoD, nD) 
{
    trackRatioSum = sum(nS)/sum(nI)
    pooledAge = (1/LamdaD) * log(LamdaD * trackRatioSum * RhoD * 
        Zeta * 0.5 + 1)/1e+06
    pooledAgeError = pooledAge * sqrt(1/sum(nI) + 1/sum(nS) + 
        1/nD + (ZetaStErr/Zeta)^2)
    output = list(pooledAge = pooledAge, StErrPooledAge = pooledAgeError)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
