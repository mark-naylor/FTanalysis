\name{BINOMFIT_PooledAge}
\alias{BINOMFIT_PooledAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the pooled age
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BINOMFIT_PooledAge(FTdataset, details = NULL, K = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FTdataset}{
%%     ~~Describe \code{FTdataset} here~~
}
  \item{details}{
%%     ~~Describe \code{details} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
}
\details{
  Used when we can assume that all of the grain have a common parent age
  --Function copied from spreadsheet calculations--
}
\value{
Returns Pooled Age +/1 1sd 
%%  If it is a LIST, use
  \item{results$age }{ Pooled age}
  \item{results$error }{ 1sd}
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
function (FTdataset, details = NULL, K = 1) 
{
    nS_Sum <- sum(FTdataset$nS)
    nI_Sum <- sum(FTdataset$nI)
    beta <- log(nS_Sum/nI_Sum)
    varBeta <- 1/nS_Sum + 1/nI_Sum
    age <- 1/LamdaD * log(1 + FTdataset$rhoD * FTdataset$Zeta * 
        LamdaD * 0.5 * nS_Sum/nI_Sum) * 1e-06
    seAge <- age * sqrt(1/nS_Sum + 1/nS_Sum + 1/FTdataset$nD + 
        (FTdataset$relErrZeta/FTdataset$Zeta)^2)
    results <- list(age = age, error = seAge)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
