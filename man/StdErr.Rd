\name{StdErr}
\alias{StdErr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
StdErr(FTdataset, PkTheta, PkFrac, PkNum, GF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{PkTheta}{
%%     ~~Describe \code{PkTheta} here~~
}
  \item{PkFrac}{
%%     ~~Describe \code{PkFrac} here~~
}
  \item{PkNum}{
%%     ~~Describe \code{PkNum} here~~
}
  \item{GF}{
%%     ~~Describe \code{GF} here~~
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
function (FTdataset, PkTheta, PkFrac, PkNum, GF) 
{
    Ns <- FTdataset$nS
    Ni <- FTdataset$nI
    Nt <- FTdataset$nT
    Num <- FTdataset$nGrain
    Zeta0 <- FTdataset$Zeta
    REZeta0 <- FTdataset$relErrZeta
    RhoD0 <- FTdataset$rhoD
    RERhoD0 <- FTdataset$relErrRhoD
    P = array(dim = c(10, 230))
    A = array(dim = c(9, 9))
    B = array(dim = c(9, 10))
    C = array(dim = c(10, 10))
    Covar = array(dim = c(19, 19))
    SEPkfrac = array(dim = (10))
    PeakAge = array(dim = (10))
    SEPeakAge = array(dim = c(10, 4))
    PkSDz = array(dim = (10))
    for (U in 1:Num) {
        SumFiu = 0
        for (I in 1:PkNum) {
            P[I, U] = PkFrac[I] * Fiu(Ns[U], Ni[U], Nt[U], PkTheta[I])
            SumFiu = SumFiu + P[I, U]
        }
        for (I in 1:PkNum) {
            P[I, U] = P[I, U]/SumFiu
        }
    }
    PIk = PkFrac[PkNum]
    for (I in 1:(PkNum - 1)) {
        for (J in 1:(PkNum - 1)) {
            A[I, J] = 0
            for (U in 1:Num) {
                Pku = P[PkNum, U]
                Term = (P[I, U]/PkFrac[I] - Pku/PIk) * (P[J, 
                  U]/PkFrac[J] - Pku/PIk)
                A[I, J] = A[I, J] + Term
            }
        }
    }
    for (I in 1:(PkNum - 1)) {
        for (J in 1:PkNum) {
            B[I, J] = 0
            for (U in 1:Num) {
                Pku = P[PkNum, U]
                if (I == J) {
                  Del1 = 1
                }
                else {
                  Del1 = 0
                }
                if (PkNum == J) {
                  Del2 = 1
                }
                else {
                  Del2 = 0
                }
                Term = P[I, U]/PkFrac[I] - Pku/PIk - Del1/PkFrac[I] + 
                  Del2/PIk
                Aju = Ns[U] - PkTheta[J] * Nt[U]
                B[I, J] = B[I, J] + P[J, U] * Aju * Term
            }
        }
    }
    for (I in 1:PkNum) {
        for (J in 1:PkNum) {
            C[I, J] = 0
            for (U in 1:Num) {
                Aiu = Ns[U] - PkTheta[I] * Nt[U]
                Aju = Ns[U] - PkTheta[J] * Nt[U]
                Term = P[J, U] * Aiu * Aju
                if (I == J) {
                  Term = Term - (Aiu * Aiu - PkTheta[I] * (1 - 
                    PkTheta[I]) * Nt[U])
                }
                C[I, J] = C[I, J] + P[I, U] * Term
            }
        }
    }
    for (I in 1:(PkNum - 1)) {
        for (J in 1:(PkNum - 1)) {
            Covar[I, J] = A[I, J]
        }
    }
    for (I in 1:(PkNum - 1)) {
        for (J in 1:PkNum) {
            Covar[I, J + PkNum - 1] = B[I, J]
        }
    }
    for (I in 1:PkNum) {
        for (J in 1:(PkNum - 1)) {
            Covar[I + PkNum - 1, J] = B[J, I]
        }
    }
    for (I in 1:PkNum) {
        for (J in 1:PkNum) {
            Covar[I + PkNum - 1, J + PkNum - 1] = C[I, J]
        }
    }
    M = 2 * PkNum - 1
    Covar[1:M, 1:M] <- solve(Covar[1:M, 1:M])
    for (I in 1:(PkNum - 1)) {
        SEPkfrac[I] = sqrt(abs(Covar[I, I]))
    }
    if (PkNum > 1) {
        Sum = 0
        for (I in 1:(PkNum - 1)) {
            for (J in 1:(PkNum - 1)) {
                print(I, J)
                Sum = Sum + Covar[I, J]
            }
        }
        SEPkfrac[PkNum] = sqrt(abs(Sum))
    }
    for (I in 1:PkNum) {
        BetaPk = log(PkTheta[I]/(1 - PkTheta[I]))
        PeakAge[I] = TaufromBeta(BetaPk, Zeta0, GF, RhoD0)
        J = I + PkNum - 1
        SE = sqrt(abs(Covar[J, J]) + RERhoD0 * RERhoD0 + REZeta0 * 
            REZeta0)
        SEPeakAge[I, 1] = TaufromBeta(BetaPk - SE, Zeta0, GF, 
            RhoD0) - PeakAge[I]
        SEPeakAge[I, 2] = TaufromBeta(BetaPk + SE, Zeta0, GF, 
            RhoD0) - PeakAge[I]
        SEPeakAge[I, 3] = TaufromBeta(BetaPk - 1.96 * SE, Zeta0, 
            GF, RhoD0) - PeakAge[I]
        SEPeakAge[I, 4] = TaufromBeta(BetaPk + 1.96 * SE, Zeta0, 
            GF, RhoD0) - PeakAge[I]
        SumNt = 0
        SumFiu = 0
        for (U in 1:Num) {
            SumNt = SumNt + P[I, U] * Nt[U]
            SumFiu = SumFiu + P[I, U]
        }
        SumNt = SumNt/SumFiu
        PkSDz[I] = sqrt(1/(SumNt * PkTheta[I] * (1 - PkTheta[I])))
    }
    Cnum <- -999
    PeakAgeCI65min = SEPeakAge[, 1]
    PeakAgeCI65plus = SEPeakAge[, 2]
    PeakAgeCI95min = SEPeakAge[, 3]
    PeakAgeCI95plus = SEPeakAge[, 4]
    results = data.frame(PeakAge, PeakAgeCI65min, PeakAgeCI65plus, 
        PeakAgeCI95min, PeakAgeCI95plus, SEPkfrac, PkSDz, Cnum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
