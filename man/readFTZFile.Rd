\name{readFTZFile}
\alias{readFTZFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readFTZFile(dir, filename, ending = ".FTZ", seperatorHeader = ",", seperatorData = "")
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
  \item{seperatorHeader}{
%%     ~~Describe \code{seperatorHeader} here~~
}
  \item{seperatorData}{
%%     ~~Describe \code{seperatorData} here~~
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
function (dir, filename, ending = ".FTZ", seperatorHeader = ",", 
    seperatorData = "") 
{
    tmp <- file(paste(dir, filename, ending, sep = ""), "rt")
    details = readLines(tmp, 1)
    close(tmp)
    header <- read.table(paste(dir, filename, ending, sep = ""), 
        skip = 1, sep = seperatorHeader, nrows = 2)
    data <- read.table(paste(dir, filename, ending, sep = ""), 
        skip = 3, sep = seperatorData)
    if (header[1, 1] == -1) {
        RERhoD = 0
        print("WARNING:: Z method!!!")
        inputFileData = NULL
        SqSize = header[1, 3]
    }
    else {
        print("Using Zeta method")
        Zeta = header[2, 1]
        StErrZeta = header[2, 2]
        REZeta = StErrZeta/Zeta
        SqSize = header[2, 3]
        RhoD = abs(header[1, 1])
        StErrRhoD = header[1, 2]
        RERhoD = StErrRhoD/100
        EffectiveUContent = header[1, 3]
        Zeta0 = Zeta
        REZeta0 = REZeta
        RhoD0 = RhoD
        RERhoD0 = RERhoD
        K = (Zeta0 * RhoD0)/(Zeta * RhoD)
        inputFileData = list(file = filename, details = details, 
            nS = data[, 1], nI = data[, 2], SqNum = data[, 3], 
            RhoD = RhoD, Zeta = Zeta, REZeta = REZeta, RERhoD = RERhoD, 
            SqSize = SqSize, StErrZeta = StErrZeta, K = K)
    }
    return(inputFileData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
