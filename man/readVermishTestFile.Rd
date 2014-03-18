\name{readVermishTestFile}
\alias{readVermishTestFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to read in the test dataset provided by Vermish with his Java RadialPlotter - for benchmarking}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readVermishTestFile(dir = "/Users/mnaylor/", filename = "FTdata_Ver.txt", ending = "", seperator = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{ending}{
%%     ~~Describe \code{ending} here~~
}
  \item{seperator}{
%%     ~~Describe \code{seperator} here~~
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
function (dir = "/Users/mnaylor/", filename = "FTdata_Ver.txt", 
    ending = "", seperator = ",") 
{
    tmp <- file(paste(dir, filename, ending, sep = ""), "rt")
    details = readLines(tmp, 1)
    close(tmp)
    header <- read.table(paste(dir, filename, ending, sep = ""), 
        skip = 1, sep = seperator, nrows = 2)
    data <- read.table(paste(dir, filename, ending, sep = ""), 
        skip = 3, sep = seperator)
    print("Using Zeta method")
    Zeta = header[1, 1]
    StErrZeta = header[1, 2]
    REZeta = StErrZeta/Zeta
    SqSize = NULL
    RhoD = abs(header[2, 1])
    StErrRhoD = header[2, 2]
    RERhoD = StErrRhoD/RhoD
    Zeta0 = Zeta
    REZeta0 = REZeta
    RhoD0 = RhoD
    RERhoD0 = RERhoD
    K = (Zeta0 * RhoD0)/(Zeta * RhoD)
    inputFileData = list(file = filename, details = details, 
        nS = data[, 1], nI = data[, 2], SqNum = data[, 3], RhoD = RhoD, 
        Zeta = Zeta, REZeta = REZeta, RERhoD = RERhoD, SqSize = SqSize, 
        StErrZeta = StErrZeta, K = K, StErrRhoD = StErrRhoD)
    return(inputFileData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
