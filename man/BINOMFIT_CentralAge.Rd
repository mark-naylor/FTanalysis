\name{BINOMFIT_CentralAge}
\alias{BINOMFIT_CentralAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard function to calculate the Central Age
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BINOMFIT_CentralAge(FTdataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
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
function (FTdataset) 
{
    nT <- FTdataset$nT
    nS <- FTdataset$nS
    nI <- FTdataset$nI
    y <- nS/nT
    z <- log10((nS + 0.5)/(nI + 0.5))
    sigma <- c()
    eta <- c()
    w <- c()
    wA <- c()
    wB <- c()
    sigma[1] <- 0.6 * sd(z)
    eta[1] <- sum(nS)/sum(nT)
    w <- nT/((eta[1] * (1 - eta[1])) + (nT - 1) * (eta[1]^2) * 
        ((1 - eta[1])^2) * sigma[1]^2)
    wA <- w * w * (y - eta[1])^2
    wB <- w * y
    for (i in 2:20) {
        eta[i] <- sum(wB)/sum(w)
        sigma[i] <- sigma[i - 1] * sqrt(sum(wA)/sum(w))
        w <- nT/((eta[i] * (1 - eta[i])) + (nT - 1) * (eta[i] * 
            eta[i]) * ((1 - eta[i]) * (1 - eta[i])) * sigma[i] * 
            sigma[i])
        wA <- w^2 * (y - eta[i])^2
        wB <- w * y
    }
    nSExp = (sum(nS) * (nS + nI))/(sum(nS) + sum(nI))
    nIExp = (sum(nI) * (nS + nI))/(sum(nS) + sum(nI))
    ChiFactor = (((nS - nSExp)^2)/nSExp) + (((nI - nIExp)^2)/nIExp)
    ChiSquared = sum(ChiFactor)
    degFreedom = length(nI) - 1
    print(degFreedom)
    centralAge <- (1/LamdaD) * (log((LamdaD * (eta[20]/(1 - eta[20])) * 
        FTdataset$rhoD * FTdataset$Zeta * 0.5) + 1))/1e+06
    centralAgeStError <- centralAge * sqrt((1/(eta[20]^2 * (1 - 
        eta[20])^2 * sum(w))) + (1/FTdataset$nD) + ((FTdataset$relErrZeta)^2))
    dispersion <- 100 * sigma[20]
    results = list(centralAge = centralAge, centralAgeStError = centralAgeStError, 
        ChiSquared = ChiSquared, degFreedom = degFreedom, eta = eta, 
        dispersion = dispersion)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
