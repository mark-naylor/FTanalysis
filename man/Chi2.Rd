\name{Chi2}
\alias{Chi2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chi2(PkTheta, PkFrac, PkNum, Num, Nt, Ns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PkTheta}{
%%     ~~Describe \code{PkTheta} here~~
}
  \item{PkFrac}{
%%     ~~Describe \code{PkFrac} here~~
}
  \item{PkNum}{
%%     ~~Describe \code{PkNum} here~~
}
  \item{Num}{
%%     ~~Describe \code{Num} here~~
}
  \item{Nt}{
%%     ~~Describe \code{Nt} here~~
}
  \item{Ns}{
%%     ~~Describe \code{Ns} here~~
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
function (PkTheta, PkFrac, PkNum, Num, Nt, Ns) 
{
    Mu = 0
    SS = 0
    for (I in 1:PkNum) {
        Mu = Mu + PkFrac[I] * PkTheta[I]
        SS = SS + PkFrac[I] * PkTheta[I] * PkTheta[I]
    }
    Factor1 = Mu * (1 - Mu)
    Factor2 = SS - Mu * Mu
    C2 = 0
    for (I in 1:Num) {
        NNt = Nt[I]
        dev = Ns[I] - Mu * NNt
        Var = Factor1 * NNt + NNt * (NNt - 1) * Factor2
        C2 = C2 + dev * dev/Var
    }
    Chi2 = C2
    Chi2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
