\name{makeFTdataset}
\alias{makeFTdataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeFTdataset(nS, nI, nD = NULL, rhoD, relErrRhoD, c, K, Zeta, relErrZeta, SqSize = NULL, geomFactor = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nS}{
%%     ~~Describe \code{nS} here~~
}
  \item{nI}{
%%     ~~Describe \code{nI} here~~
}
  \item{nD}{
%%     ~~Describe \code{nD} here~~
}
  \item{rhoD}{
%%     ~~Describe \code{rhoD} here~~
}
  \item{relErrRhoD}{
%%     ~~Describe \code{relErrRhoD} here~~
}
  \item{c}{
%%     ~~Describe \code{c} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{Zeta}{
%%     ~~Describe \code{Zeta} here~~
}
  \item{relErrZeta}{
%%     ~~Describe \code{relErrZeta} here~~
}
  \item{SqSize}{
%%     ~~Describe \code{SqSize} here~~
}
  \item{geomFactor}{
%%     ~~Describe \code{geomFactor} here~~
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
function (nS, nI, nD = NULL, rhoD, relErrRhoD, c, K, Zeta, relErrZeta, 
    SqSize = NULL, geomFactor = 0.5) 
{
    if (length(nI) == length(nS)) {
        nGrain = length(nS)
    }
    else {
        print("ERROR: There are different numbers of nI and nS")
    }
    NNi = K * nI
    nT = nS + NNi
    NNs = nS
    NNi = NNi
    b <- Zeta * rhoD * 0.5 * LamdaD
    Zgrain = log(b * NNs/NNi)
    Zerr = sqrt(1/NNs + 1/NNi)
    QMeanZerr = sqrt(sum(Zerr * Zerr)/nGrain)
    stErrRhoD <- relErrRhoD * rhoD
    nD <- (rhoD/stErrRhoD)^2
    grainAges = TaufromZ(Zgrain)
    grainAges95min = TaufromZ(Zgrain + 2 * Zerr)
    grainAges95max = TaufromZ(Zgrain - 2 * Zerr)
    FTdataset = list(nS = nS, nI = nI, nD = nD, nT = nT, grainAges = grainAges, 
        rhoD = rhoD, relErrRhoD = relErrRhoD, Zeta = Zeta, relErrZeta = relErrZeta, 
        c = c, nGrain = nGrain, Zgrain = Zgrain, Zerr = Zerr, 
        QMeanZerr = QMeanZerr, b = b, nD = nD, grainAges95min = grainAges95min, 
        grainAges95max = grainAges95max)
    class(FTdataset) <- "FTdataset"
    return(FTdataset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
