\name{TaufromBeta}
\alias{TaufromBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TaufromBeta(Beta, Zeta0, GF, RhoD0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Beta}{
%%     ~~Describe \code{Beta} here~~
}
  \item{Zeta0}{
%%     ~~Describe \code{Zeta0} here~~
}
  \item{GF}{
%%     ~~Describe \code{GF} here~~
}
  \item{RhoD0}{
%%     ~~Describe \code{RhoD0} here~~
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
function (Beta, Zeta0, GF, RhoD0) 
{
    if (Beta > 80) {
        BBeta <- 80
    }
    else if (Beta < -80) {
        BBeta <- -80
    }
    else {
        BBeta <- Beta
    }
    ExpZm = LamdaD * Zeta0 * GF * RhoD0 * exp(BBeta)
    TaufromBeta = log(1 + ExpZm)/(1e+06 * LamdaD)
    TaufromBeta
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
