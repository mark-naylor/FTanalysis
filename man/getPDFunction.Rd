\name{getPDFunction}
\alias{getPDFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for calculating the PD function for the PD plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPDFunction(FTdataset, zeroNsOffset = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{zeroNsOffset}{
%%     ~~Describe \code{zeroNsOffset} here~~
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
function (FTdataset, zeroNsOffset = 0) 
{
    if (zeroNsOffset == 0) {
        Zgrain = FTdataset$Zgrain
        Zerr = FTdataset$Zerr
    }
    else {
        b <- FTdataset$b
        NNs <- FTdataset$nS + zeroNsOffset
        NNi <- FTdataset$nI + zeroNsOffset
        Zgrain = log(b * NNs/NNi)
        Zerr = sqrt(1/NNs + 1/NNi)
    }
    Zmax = max(Zgrain)
    Zmin = min(Zgrain)
    Num = length(Zgrain)
    QMeanZerr = sum(Zerr * Zerr)
    QMeanZerr = sqrt(QMeanZerr/(Num))
    barWidth = 0.1
    zWidth = barWidth/5
    LwLmt = barWidth * (floor(Zmin/barWidth) - 1)
    UpLmt = barWidth * (1 + floor(Zmax/barWidth))
    LwLmt = Zmin - (5 * QMeanZerr)
    UpLmt = Zmax + (5 * QMeanZerr)
    Zi = LwLmt
    atZ = c()
    gCount = floor((UpLmt - LwLmt)/zWidth) + 1
    pd = array(dim = c(gCount, 4))
    for (j in 1:gCount) {
        tmp <- KernelPD(Zi, -1, Zgrain, Zerr, Num)
        Dsum <- tmp[1]
        SEDsum <- tmp[2]
        Dsum = 100 * barWidth * Dsum/Num
        SEDsum = 100 * barWidth * SEDsum/Num
        pd[j, 1] = TaufromZ(Zi)
        pd[j, 2] = Dsum
        pd[j, 3] = Dsum - SEDsum
        pd[j, 4] = Dsum + SEDsum
        atZ[j] = Zi
        Zi = Zi + zWidth
    }
    dat2 = data.frame(atZ, pd[, 1], pd[, 2], pd[, 3], pd[, 4], 
        log(pd[, 1]))
    colnames(dat2) <- c("Zi", "age", "mean", "upper", "lower", 
        "logAges")
    return(dat2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
