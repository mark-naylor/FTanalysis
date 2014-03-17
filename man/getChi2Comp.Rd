\name{getChi2Comp}
\alias{getChi2Comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getChi2Comp(Chi1, DF1, Chi2, DF2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Chi1}{
%%     ~~Describe \code{Chi1} here~~
}
  \item{DF1}{
%%     ~~Describe \code{DF1} here~~
}
  \item{Chi2}{
%%     ~~Describe \code{Chi2} here~~
}
  \item{DF2}{
%%     ~~Describe \code{DF2} here~~
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
function (Chi1, DF1, Chi2, DF2) 
{
    print("===================Chi2Comp Program v 1.1  (Brandon 11/12/93)===================")
    print("Compares X^2 values to determine if fit is significantly")
    print("improved. See Bevington (1969, p. 200) for description of method.")
    print("The routine for calculating the F distribution probability is from Press")
    print("and others (1986).")
    print("  ")
    print("F = [Chi^2(1) - Chi^2(2)] / [Chi^2(2)/DF2]")
    print("    where 1 refers to the first fit and 2,")
    print("    the second fit with Chi^2(1) >= Chi^2(2).")
    print("  ")
    print("The 'second fit becomes first fit' option allows the user to efficiently")
    print("analyze a series of results where X^2 is successively reduced due to")
    print("an increase in the number of fit parameters. This option is particularly")
    print("useful for examining data from the BINOMFIT AND GAUSSFIT programs.")
    print("To halt program, use Ctrl-Break.")
    print(paste("First: Chi^2, degrees of freedom: ", Chi1, DF1))
    print(paste("Second: Chi^2, degrees of freedom: ", Chi2, 
        DF2))
    F = (Chi1 - Chi2)/(Chi2/DF2)
    print(paste("F = ", F))
    DF1 = DF1 - DF2
    X = DF2/(DF2 + DF1 * F)
    Prob = Ibeta(a = DF2/2, b = DF1/2, x = X)
    print(paste("Probability that F is by chance alone (\%) ", 
        100 * Prob))
    Prob
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
