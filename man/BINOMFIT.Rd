\name{BINOMFIT}
\alias{BINOMFIT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BINOMFIT(FTdataset, TrialAges, PkNum, details = NULL, K = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{TrialAges}{
%%     ~~Describe \code{TrialAges} here~~
}
  \item{PkNum}{
%%     ~~Describe \code{PkNum} here~~
}
  \item{details}{
%%     ~~Describe \code{details} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (FTdataset, TrialAges, PkNum, details = NULL, K = 1) 
{
    PeakAge = array(dim = (10))
    SEPeakAge = array(dim = c(10, 4))
    PkTheta = array(dim = (10))
    PkFrac = array(dim = (10))
    SEPkfrac = array(dim = (10))
    PkSDz = array(dim = (10))
    RhoD = FTdataset$rhoD
    relErrRhoD = FTdataset$relErrRhoD
    Zeta = FTdataset$Zeta
    relErrZeta = FTdataset$relErrZeta
    Zgrain = FTdataset$Zgrain
    Zerr = FTdataset$Zerr
    Ns = FTdataset$nS
    Ni = FTdataset$nI
    Nt = FTdataset$nT
    grainAges = FTdataset$grainAges
    QMeanZerr = FTdataset$QMeanZerr
    nGrain = FTdataset$nGrain
    peak = array(dim = c(PkNum, 3))
    GF = 0.5
    KernelFactor = 0.6
    PkWidthFactor = sqrt(1 + KernelFactor * KernelFactor)
    Label1 = "=================BinomFit Program  v. 1.8  (Brandon 9/15/04)=================="
    print(Label1)
    print("This program determines best-fit peaks based on a binomial model described by")
    print("Galbraith and Green (1990, Nuclear Tracks and Radiation Measurements, v. 17,")
    print("p. 197-206). The user must give an initial guess for the number of peaks")
    print("and the approximate ages of the peaks. The fit is found using the EM method")
    print("outlined in Galbraith (1988, Technometrics, v. 30, p. 280).")
    print("At present, the program can fit up to 10 peaks.")
    print(" ")
    print("Note that BINOMFIT does not include the routine used in ZETAAGE to ")
    print("calculate  precise ages for grains with low track densities. As a result,")
    print("the range in grain ages reported here may differ somewhat from that reported")
    print("in ZETAAGE but this difference has no influence on the BINOMFIT results.")
    print(" ")
    print("The input data should use the same format as for the ZETAAGE program:")
    print("Line 1:  Title for data set")
    print("Line 2:  This line is flagged by making the first value negative.")
    print("         All values on this line should be separated by commas. ")
    print("     ZETA METHOD: Effective track density (tr/cm^2), the relative")
    print("         standard error (\%) for the fluence monitor, and the effective ")
    print("         uranium concentration (ppm) of the monitor standard. ")
    print("     Z METHOD: Enter line as:  -1, 0, 0.")
    print("Line 3:  All values should be separated by commas. ")
    print("     ZETA METHOD: Zeta factor (yr cm^2/tr), standard error (yr cm^2/tr), ")
    print("         counter square size (cm^2).")
    print("     Z METHOD: Z factor (yr), standard error (yr), counter square size (cm^2).")
    print("Line 4 and on: For each grain, enter the measured spontaneous and induced ")
    print("         tracks, and the number of squares for the area counted. Data ")
    print("         values can be separated by spaces or commas.")
    print("NOTES: To merge another dataset, repeat lines 2 and on.")
    print("Input routine skips comment lines (first character=semicolon) and blank lines.")
    print("Methods and fission track constants are after Hurford and Green (1983).")
    print(paste("Total decay constant for 238U (yr^-1) = ", LamdaD))
    print(paste("Geometry factor = ", GF))
    print("")
    Zmax = max(Zgrain)
    Zmin = min(Zgrain)
    grainAgeMin = min(grainAges)
    grainAgeMax = max(grainAges)
    print(grainAges)
    EstPkWidth = PkWidthFactor * QMeanZerr
    PeakAge[1:PkNum] = TrialAges[1:PkNum]
    PkTheta = Theta(PeakAge, Zeta, GF, RhoD)
    for (I in 1:PkNum) {
        Zpk = ZfromTau(PeakAge[I])
        PkDen = SumPD(Zpk, 0, Zgrain, Zerr, nGrain, KernelFactor)
        PkFrac[I] = PkDen * 0.15/(C0 * nGrain)
        if (PkFrac[I] < 0.001) {
            PkFrac[I] = 0.001
        }
    }
    print(" ")
    print(Label1)
    print("FILE HEADER:")
    print(details)
    print("")
    print("FIT OPTION: Best-fit peaks using the binomial model of Galbraith and Green")
    print("")
    print("---------------------INITIAL GUESS for MODEL PARAMETERS------------------------")
    print(paste("NUMBER OF PEAKS TO FIT = ", PkNum))
    print("PEAK #)  PEAK AGE     THETA    FRACTION(\%)  COUNT")
    for (I in 1:PkNum) {
        print(paste(I, PeakAge[I], PkTheta[I], PkFrac[I] * 100, 
            nGrain * PkFrac[I]))
    }
    print(paste("                    TOTAL RANGE for GRAIN AGES = ", 
        grainAgeMin, " to ", grainAgeMax, " Ma"))
    print(paste("                              NUMBER OF GRAINS = ", 
        nGrain))
    print(paste("                    DEGREES OF FREEDOM for FIT = ", 
        nGrain - (2 * PkNum - 1)))
    print(paste("         AVERAGE OF THE SE(Z)'S for THE GRAINS = ", 
        QMeanZerr))
    print(paste("ESTIMATED WIDTH OF PEAKS IN PD PLOT IN Z UNITS = ", 
        EstPkWidth))
    print("   *** EM routine is finding a set of best-fit peaks ***")
    results <- EM(PkTheta, PkFrac, PkNum, nGrain, Ns, Ni, Nt)
    results2 <- StdErr(FTdataset = FTdataset, PkTheta = results$PkTheta, 
        PkFrac = results$PkFrac, PkNum = PkNum, GF = GF)
    print("")
    print("------------------------PARAMETERS for BEST-FIT PEAKS--------------------------")
    print("{Standard error for peak age includes group error}")
    print(paste("{Peak width is for PD plot assuming a kernel factor = ", 
        KernelFactor, "}"))
    print("PEAK   ----PEAK AGE & CONF. INTERVAL (MA)---- PEAK WIDTH --GRAINS IN PEAK----")
    print("NUMBER   MEAN (----68\% CI---) (----95\% CI---)    W(Z)    FRAC(\%) SE(\%)  COUNT")
    TotalNum = 0
    TotalFrac = 0
    D1 <- array(dim = PkNum)
    D2 <- array(dim = PkNum)
    D3 <- array(dim = PkNum)
    for (I in 1:PkNum) {
        TotalNum = TotalNum + results$PkFrac[I] * nGrain
        TotalFrac = TotalFrac + results$PkFrac[I]
        D1[I] = 100 * results$PkFrac[I]
        D2[I] = 100 * results2$SEPkfrac[I]
        D3[I] = results$PkFrac[I] * nGrain
        PkSDz[I] = PkWidthFactor * results2$PkSDz[I]
        print(paste(I, results2$PeakAge[I], results2$PeakAgeCI65min[I], 
            results2$PeakAgeCI65plus[I], results2$PeakAgeCI95min[I], 
            results2$PeakAgeCI95plus[I], PkSDz[I], D1[I], D2[I], 
            D3[I], sep = "   "))
    }
    print(paste("Total:", 100 * TotalFrac, TotalNum, sep = "   "))
    print(" ")
    ChiSq = Chi2(results$PkTheta, results$PkFrac, PkNum, nGrain, 
        Nt, Ns)
    RChiSq = ChiSq/(nGrain - (2 * PkNum - 1))
    print(paste("        LOG-LIKELIHOOD for BEST FIT = ", results$LogLike[1]))
    print(paste("     CHI-SQUARED VALUE for BEST FIT = ", ChiSq))
    print(paste("          REDUCED CHI-SQUARED VALUE = ", RChiSq))
    print(paste("         DEGREES OF FREEDOM for FIT = ", nGrain - 
        (2 * PkNum - 1)))
    print(paste("               NUMBER OF ITERATIONS = ", results$Iter[1]))
    resultsOutput <- list(PeakAgeResults = results2$PeakAge, 
        PeakAgeCI65min = results2$PeakAgeCI65min, PeakAgeCI65plus = results2$PeakAgeCI65plus, 
        PeakAgeCI95min = results2$PeakAgeCI95min, PeakAgeCI95plus = results2$PeakAgeCI95plus, 
        PkSDz = PkSDz, D1 = D1, D2 = D2, D3 = D3, PeakAgeProposed = PeakAge, 
        PkSDz = PkSDz, PkFrac = results$PkFrac, PkNum = PkNum, 
        PkZ = results$PkZ, logLike = results$LogLike[1], ChiSq = ChiSq, 
        ReducedChiSq = RChiSq, degFreedom = (nGrain - (2 * PkNum - 
            1)))
    return(resultsOutput)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
