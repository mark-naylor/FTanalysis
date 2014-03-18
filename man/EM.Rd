\name{EM}
\alias{EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expectation Maximisation algorithm used for the binomfit peak fitting 
}
\description{
(modified from Binomfit)
}
\usage{
EM(PkTheta, PkFrac, PkNum, Num, Ns, Ni, Nt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PkTheta}{
%%     ~~Describe \code{PkTheta} here~~
}
  \item{PkFrac}{
%%     ~~Describe \code{PkFrac} here~~
}
  \item{PkNum}{
%%     ~~Describe \code{PkNum} here~~
}
  \item{Num}{
%%     ~~Describe \code{Num} here~~
}
  \item{Ns}{
%%     ~~Describe \code{Ns} here~~
}
  \item{Ni}{
%%     ~~Describe \code{Ni} here~~
}
  \item{Nt}{
%%     ~~Describe \code{Nt} here~~
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
function (PkTheta, PkFrac, PkNum, Num, Ns, Ni, Nt) 
{
    Rtol = 1e-05
    P = array(dim = c(10, 230))
    print("ITERATION  LOG-LIKELIHOOD  RELATIVE CHANGE")
    OldLogLike = 0
    for (Iter in 1:300) {
        LogLike = 0
        for (U in 1:Num) {
            SumFiu = 0
            for (I in 1:PkNum) {
                P[I, U] = PkFrac[I] * Fiu(Ns[U], Ni[U], Nt[U], 
                  PkTheta[I])
                SumFiu = SumFiu + P[I, U]
            }
            for (I in 1:PkNum) {
                P[I, U] = P[I, U]/SumFiu
            }
            LogLike = LogLike + log(SumFiu)
        }
        for (I in 1:PkNum) {
            SumPkfrac = 0
            for (U in 1:Num) {
                SumPkfrac = SumPkfrac + P[I, U]
            }
            PkFrac[I] = SumPkfrac/Num
        }
        for (I in 1:PkNum) {
            Sum1 = 0
            Sum2 = 0
            for (U in 1:Num) {
                Sum1 = Sum1 + P[I, U] * Ns[U]
                Sum2 = Sum2 + P[I, U] * Nt[U]
            }
            PkTheta[I] = Sum1/Sum2
        }
        ChiSq = Chi2(PkTheta, PkFrac, PkNum, Num, Nt, Ns)
        Rchange = abs(0.5 * (LogLike - OldLogLike)/(LogLike + 
            OldLogLike))
        print(paste(Iter, LogLike, Rchange))
        if (Rchange <= Rtol) {
            Itest1 = Itest1 + 1
            if (Itest1 == 3) {
                break
            }
        }
        else {
            Itest1 = 0
        }
        OldLogLike = LogLike
        if (Iter > 300) {
            print()
            print("   >>>    WARNING: LOCAL MINIMUM POORLY DEFINED   <<<")
            print("   >>> CHECK THAT BEST-FIT SOLUTION IS REASONABLE <<<")
            break
        }
    }
    results = data.frame(PkTheta, PkFrac, LogLike, Iter)
    results
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
