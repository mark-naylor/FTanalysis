\name{Fiu}
\alias{Fiu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the Probability Mass Function for use in the EM algorithm for estimating peak ages.
}
\description{
F_iu(y_u;beta_i) is the probability of the observation y_u give the age model beta_i
}
\usage{
Fiu(Ru, Su, Tu, PkTheta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ru}{
%%     ~~Describe \code{Ru} here~~
}
  \item{Su}{
%%     ~~Describe \code{Su} here~~
}
  \item{Tu}{
%%     ~~Describe \code{Tu} here~~
}
  \item{PkTheta}{
%%     ~~Describe \code{PkTheta} here~~
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
function (Ru, Su, Tu, PkTheta) 
{
    if ((Tu > 100) && (PkTheta > 0.05) && (PkTheta < 0.95)) {
        Mean = Tu * PkTheta
        SD = sqrt(Mean * (1 - PkTheta))
        Z = (Ru - Mean)/SD
        Tmp = (C0/SD) * exp(-0.5 * Z * Z)
    }
    else {
        Tmp = Bico(Tu, Ru) * (PkTheta^Ru) * (1 - PkTheta)^Su
    }
    if (Tmp <= 1e-10) {
        Fiu = 1e-10
    }
    else {
        Fiu = Tmp
    }
    Fiu
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
