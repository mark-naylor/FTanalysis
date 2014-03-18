\name{Ftest}
\alias{Ftest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform the Ftest comparing the Chi2 and df from 2 different age models on the same data
}
\description{
 ===================Chi2Comp Program v 1.1  (Brandon 11/12/93)===================
Compares X^2 values to determine if fit is significantly
improved. See Bevington (1969, p. 200) for description of method.

}
\usage{
Ftest(Chi1, Chi2, DF1, DF2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Chi1}{
%%     ~~Describe \code{Chi1} here~~
}
  \item{Chi2}{
Chi^2(1) >= Chi^2(2)
  \item{DF1}{
%%     ~~Describe \code{DF1} here~~
}
  \item{DF2}{
%%     ~~Describe \code{DF2} here~~
}
}
\details{
The routine for calculating the F distribution probability is from Press
and others (1986)

F = [Chi^2(1) - Chi^2(2)] / [Chi^2(2)/DF2]
    where 1 refers to the first fit and 2
    the second fit with Chi^2(1) >= Chi^2(2)
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
function (Chi1, Chi2, DF1, DF2) 
{
    if (Chi2 > Chi1) {
        tmp1 <- Chi1
        tmp2 <- DF1
        Chi1 <- Chi2
        DF1 <- DF2
        Chi2 <- tmp1
        DF2 <- tmp2
    }
    F <- (Chi1 - Chi2)/(Chi2/DF2)
    DFdiff <- DF1 - DF2
    X <- DF2/(DF2 + DFdiff * F)
    Prob <- Ibeta(a = DF2/2, b = DFdiff/2, x = X)
    return(Prob)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
