\name{KernelPD}
\alias{KernelPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CHECK binomfit description
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KernelPD(Zi, Order, Zgrain, Zerr, Num, KernelFactor = 0.6)
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
    PD = 0
    SSPD = 0
    for (I in 1:Num) {
        Dev = Zi - Zgrain[I]
        s = KernelFactor * Zerr[I]
        if (abs(Dev) < (5 * s)) {
            Dev2 = Dev * Dev
            S2 = s * s
            ExpTerm = exp(-Dev2/(2 * S2))
            if (Order == -1) {
                Term = (C0/s) * ExpTerm
                PD = PD + Term
                SSPD = SSPD + Term * Term
            }
            else if (Order == 0) {
                PD = PD + (C0/s) * ExpTerm
            }
            else if (Order == 1) {
                S3 = s * S2
                PD = PD - (C0 * Dev/S3) * ExpTerm
            }
            else if (Order == 2) {
                S5 = S2 * S2 * s
                PD = PD + (C0 * (Dev2 - S2)/S5) * ExpTerm
            }
            else if (Order == 3) {
                S7 = S2 * S2 * S2 * s
                PD = PD + (C0 * Dev * (-Dev2 + 3 * S2)/S7) * 
                  ExpTerm
            }
        }
    }
    if ((Order == -1) && (Num > 1)) {
        SEPD = sqrt((Num * SSPD - PD * PD)/(Num - 1))
    }
    else {
        SEPD = 0
    }
    results <- c(PD, SEPD)
    results
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
