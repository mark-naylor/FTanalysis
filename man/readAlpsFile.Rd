\name{readAlpsFile}
\alias{readAlpsFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readAlpsFile(i = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
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
function (i = 2) 
{
    dataDir = "/Users/mnaylor/Documents/MyPapers/PaperSubmitted/SinclairAndNaylor_2010_bias/AlpineRawData/"
    dataFiles = c("99-11", "C7-2(1)", "E2-1A(1)", "N13-1(1)", 
        "N6-2(1)", "99-12", "C9-2(1)", "E2-1C(1)", "N13-3(1)", 
        "99-22", "E1-2(1)", "E2-2A(1)", "N5-2(1)", "C10-2(1)", 
        "E18-2(1)", "N12-1B(1)", "N5-4(1)")
    file = dataFiles[i]
    inputData = readBinimfitOutputAsCSV(dataDir, file)
    return(inputData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
