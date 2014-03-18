\name{Theta}
\alias{Theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Given a FT age, this routine calculates Theta = Ns /(Ns+Ni) each grain to be used in the EM algorithm}
\description{
}
\usage{
Theta(age, Zeta0, GF, RhoD0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{
FT age in Ma
}
  \item{Zeta0}{
zeta factor
}
  \item{GF}{
Geometry factor
}
  \item{RhoD0}{
monitor track density used for this dataset
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
function (age, Zeta0, GF, RhoD0) 
{
    ExpTL = exp(age * 1e+06 * LamdaD)
    Theta = (ExpTL - 1)/(LamdaD * Zeta0 * GF * RhoD0 + ExpTL - 
        1)
    Theta
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
