\name{PDplot}
\alias{PDplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PDplot(FTdataset, resultsOutput, plotType = 1, zeroNsOffset = 0, ageLabels = c(1, 2, 3, 5, 7, 10, 20, 50, 100, 200))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{resultsOutput}{
%%     ~~Describe \code{resultsOutput} here~~
}
  \item{plotType}{
%%     ~~Describe \code{plotType} here~~
}
  \item{zeroNsOffset}{
%%     ~~Describe \code{zeroNsOffset} here~~
}
  \item{ageLabels}{
%%     ~~Describe \code{ageLabels} here~~
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
function (FTdataset, resultsOutput, plotType = 1, zeroNsOffset = 0, 
    ageLabels = c(1, 2, 3, 5, 7, 10, 20, 50, 100, 200)) 
{
    PeakAgeResults = resultsOutput$PeakAgeResults
    PeakAgeProposed = resultsOutput$PeakAgeProposed
    PkSDz = resultsOutput$PkSDz
    PkFrac = resultsOutput$PkFrac
    PkNum = resultsOutput$PkNum
    PkZ = resultsOutput$PkZ
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
    PkZ = ZfromTau(PeakAgeResults)
    Zmax = max(Zgrain)
    Zmin = min(Zgrain)
    Ages = TaufromZ(Zgrain)
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
    pd = array(dim = c(gCount, 5 + PkNum))
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
        SumPkD = 0
        for (i in 1:PkNum) {
            dev = Zi - PkZ[i]
            S = PkSDz[i]
            if (abs(dev) < 5 * S) {
                PkD = (100 * barWidth * PkFrac[i]) * (C0/S) * 
                  exp(-(dev * dev)/(2 * S * S))
            }
            else {
                PkD = 0
            }
            SumPkD = SumPkD + PkD
            pd[j, 4 + i] = PkD
        }
        pd[j, 5 + i] = SumPkD
        Zi = Zi + zWidth
    }
    dat = data.frame(Ages, Zgrain, log(Ages))
    colnames(dat) <- c("Ages", "Zgrain", "logAges")
    dat2 = data.frame(atZ, pd[, 1], pd[, 2], pd[, 3], pd[, 4], 
        log(pd[, 1]))
    colnames(dat2) <- c("Zi", "age", "mean", "upper", "lower", 
        "logAges")
    bottom = pd[, 1] * 0
    if (PkNum == 1) {
        dat3 = data.frame(atZ, pd[, 1], pd[, 5], pd[, 6], log(pd[, 
            1]), bottom)
        colnames(dat3) <- c("Zi", "age", "firstPeak", "total", 
            "logAges", "bottom")
    }
    else if (PkNum == 2) {
        dat3 = data.frame(atZ, pd[, 1], pd[, 5], pd[, 6], pd[, 
            7], log(pd[, 1]), bottom)
        colnames(dat3) <- c("Zi", "age", "firstPeak", "secondPeak", 
            "total", "logAges", "bottom")
    }
    else if (PkNum == 3) {
        dat3 = data.frame(atZ, pd[, 1], pd[, 5], pd[, 6], pd[, 
            7], pd[, 8], log(pd[, 1]), bottom)
        colnames(dat3) <- c("Zi", "age", "firstPeak", "secondPeak", 
            "thirdPeak", "total", "logAges", "bottom")
    }
    else if (PkNum == 4) {
        dat3 = data.frame(atZ, pd[, 1], pd[, 5], pd[, 6], pd[, 
            7], pd[, 8], pd[, 9], log(pd[, 1]), bottom)
        colnames(dat3) <- c("Zi", "age", "firstPeak", "secondPeak", 
            "thirdPeak", "fourthPeak", "total", "logAges", "bottom")
    }
    if (plotType == 1) {
        p <- ggplot(dat, aes(Ages)) + geom_histogram(binwidth = barWidth) + 
            scale_x_continuous(trans = "log", breaks = ageLabels, 
                expand = c(0.001, 0))
    }
    else if (plotType == 2) {
        p <- ggplot(data = dat2) + geom_line(aes(x = age, y = mean)) + 
            geom_line(aes(x = age, y = lower), colour = "red", 
                linetype = "dashed") + geom_line(aes(x = age, 
            y = upper), linetype = "dashed", colour = "red") + 
            xlim(0, 50)
    }
    else if (plotType == 3) {
        p <- ggplot(dat, aes(Ages)) + geom_histogram(binwidth = 1) + 
            geom_line(aes(x = age, y = firstPeak), data = dat3, 
                colour = "red", linetype = "dotdash") + geom_line(aes(x = age, 
            y = secondPeak), data = dat3, colour = "purple", 
            linetype = "dotdash") + geom_line(aes(x = age, y = thirdPeak), 
            data = dat3, colour = "yellow", linetype = "dotdash") + 
            geom_line(aes(x = age, y = total), data = dat3, colour = "red") + 
            xlim(0, 50)
    }
    else if (plotType == 4) {
        p <- ggplot(data = dat3) + geom_line(aes(x = age, y = mean), 
            colour = "red", linetype = "dotdash") + geom_line(aes(x = age, 
            y = lower), colour = "purple", linetype = "dotdash") + 
            geom_line(aes(x = age, y = upper), colour = "blue", 
                linetype = "dotdash") + geom_line(aes(x = age, 
            y = total), colour = "black") + xlim(0, 50)
    }
    else if (plotType == 5) {
        p <- ggplot(dat, aes(Zgrain)) + geom_histogram(binwidth = barWidth) + 
            geom_line(aes(x = Zi, y = mean), data = dat2, colour = "red") + 
            geom_line(aes(x = Zi, y = lower), data = dat2, colour = "red", 
                linetype = "dashed") + geom_line(aes(x = Zi, 
            y = upper), data = dat2, linetype = "dashed", colour = "red")
    }
    else if (plotType == 6) {
        lims = ageLabels
        dfLims = data.frame(x = lims)
        p <- ggplot(dat, aes(Ages)) + geom_histogram(binwidth = barWidth) + 
            scale_x_continuous(trans = "log", breaks = ageLabels, 
                expand = c(0.001, 0)) + geom_line(aes(x = age, 
            y = mean), data = dat2, colour = "red") + geom_line(aes(x = age, 
            y = lower), data = dat2, colour = "red", linetype = "dashed") + 
            geom_line(aes(x = age, y = upper), data = dat2, linetype = "dashed", 
                colour = "red") + xlab("GrainAge [Ma]") + ylab("Count") + 
            ggtitle("PD plot with grain age density")
    }
    else if (plotType == 7) {
        lims = c(1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 120, 
            130, 150, 170, 200)
        dfLims = data.frame(x = lims)
        p <- ggplot(dat, aes(Ages)) + geom_histogram(binwidth = barWidth) + 
            scale_x_continuous(trans = "log", breaks = ageLabels, 
                expand = c(0.001, 0)) + geom_ribbon(aes(x = age, 
            ymax = firstPeak, ymin = bottom), data = dat3, fill = "purple", 
            alpha = 0.3, colour = "purple", linetype = "dotdash", 
            size = 0.8) + geom_line(aes(x = age, y = firstPeak), 
            data = dat3, colour = "purple", linetype = "dotdash", 
            size = 0.8)
        if (PkNum > 1) {
            p <- p + geom_ribbon(aes(x = age, ymax = secondPeak, 
                ymin = bottom), data = dat3, fill = "green", 
                alpha = 0.3, colour = "green", linetype = "dotdash", 
                size = 0.8) + geom_line(aes(x = age, y = secondPeak), 
                data = dat3, colour = "green", linetype = "dotdash", 
                size = 0.8)
        }
        if (PkNum > 2) {
            p <- p + geom_ribbon(aes(x = age, ymax = thirdPeak, 
                ymin = bottom), data = dat3, fill = "yellow", 
                alpha = 0.3, colour = "yellow", linetype = "dotdash", 
                size = 0.8) + geom_line(aes(x = age, y = thirdPeak), 
                data = dat3, colour = "yellow", linetype = "dotdash", 
                size = 0.8)
        }
        if (PkNum > 3) {
            p <- p + geom_ribbon(aes(x = age, ymax = fourthPeak, 
                ymin = bottom), data = dat3, fill = "red", alpha = 0.3, 
                colour = "red", linetype = "dotdash", size = 0.8) + 
                geom_line(aes(x = age, y = fourthPeak), data = dat3, 
                  colour = "red", linetype = "dotdash", size = 0.8)
        }
        p <- p + geom_line(aes(x = age, y = total), data = dat3, 
            colour = "red", size = 0.9) + xlab("GrainAge [Ma]") + 
            ylab("Count") + ggtitle("PD plot with inverted age components")
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
