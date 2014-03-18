\name{SumPD}
\alias{SumPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SumPD(Zi, Order, Zgrain, Zerr, Num, KernelFactor = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Zi}{
%%     ~~Describe \code{Zi} here~~
}
  \item{Order}{
%%     ~~Describe \code{Order} here~~
}
  \item{Zgrain}{
%%     ~~Describe \code{Zgrain} here~~
}
  \item{Zerr}{
%%     ~~Describe \code{Zerr} here~~
}
  \item{Num}{
%%     ~~Describe \code{Num} here~~
}
  \item{KernelFactor}{
%%     ~~Describe \code{KernelFactor} here~~
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
function (Zi, Order, Zgrain, Zerr, Num, KernelFactor = 0.6) 
{
    Dsum = 0
    for (I in 1:Num) {
        dev = Zi - Zgrain[I]
        S = KernelFactor * Zerr[I]
        if (abs(dev) < 5 * S) {
            Dev2 = dev * dev
            S2 = S * S
            ExpTerm = exp(-Dev2/(2 * S2))
            if (Order == 0) {
                Dsum = Dsum + (C0/S) * ExpTerm
            }
            else if (Order == 1) {
                S3 = S * S2
                Dsum = Dsum - (C0 * dev/S3) * ExpTerm
            }
            else if (Order == 2) {
                S5 = S2 * S2 * S
                Dsum = Dsum + (C0 * (Dev2 - S2)/S5) * ExpTerm
            }
            else if (Order == 3) {
                S7 = S2 * S2 * S2 * S
                Dsum = Dsum + (C0 * dev * (-Dev2 + 3 * S2)/S7) * 
                  ExpTerm
            }
        }
    }
    SumPD = Dsum
    SumPD
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
