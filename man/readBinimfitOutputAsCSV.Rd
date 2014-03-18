\name{readBinimfitOutputAsCSV}
\alias{readBinimfitOutputAsCSV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to read in the standard binomfit output files to enable analysis of existing datasets and back catalogues
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readBinimfitOutputAsCSV(dataDir, dataFilename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataDir}{
%%     ~~Describe \code{dataDir} here~~
}
  \item{dataFilename}{
%%     ~~Describe \code{dataFilename} here~~
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
function (dataDir, dataFilename) 
{
    filename = paste(dataDir, dataFilename, ".csv", sep = "")
    paras = scan(filename, skip = 10, nlines = 1, what = list(x = "", 
        y = 0, z = ""), sep = ",")
    ZetaType = paras$x[1]
    Zeta = paras$y[1]
    paras = scan(filename, skip = 11, nlines = 1, what = list(x = "", 
        y = 0, z = "", z = "", z = "", z = "", z = "", z = "", 
        z = ""), sep = ",")
    StErrZeta = paras$y[1]
    REZeta = StErrZeta/Zeta
    paras = scan(filename, skip = 12, nlines = 1, what = list(x = "", 
        y = 0, z = "", z = "", z = "", z = "", z = "", z = "", 
        z = ""), sep = ",")
    nGrains = paras$y[1]
    paras = scan(filename, skip = 16, nlines = 1, what = list(x = "", 
        y = 0, z = "", z = "", z = "", z = "", z = "", z = "", 
        z = ""), sep = ",")
    Nd = paras$y[1]
    paras = scan(filename, skip = 17, nlines = 1, what = list(x = "", 
        y = 0, z = "", z = "", z = "", z = "", z = "", z = "", 
        z = ""), sep = ",")
    RhoD = paras$y[1]
    RERhoD = 0
    K = 1
    data = read.csv(file = filename, skip = 22, header = TRUE, 
        sep = ",", nrows = nGrains)
    nI = data$Ni
    nS = data$Ns
    grainAges = data$age
    inputFileData = list(nI = nI, nS = nS, grainAges = grainAges, 
        RhoD = RhoD, RERhoD = RERhoD, Zeta = Zeta, REZeta = REZeta, 
        StErrZeta = StErrZeta, Nd = Nd, K = K)
    return(inputFileData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
