\name{plot4x2AgeHists}
\alias{plot4x2AgeHists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot4x2AgeHists(ageList1, ageList2, ageList3, ageList5, preferedAgeModels1, preferedAgeModels2, preferedAgeModels3, preferedAgeModels4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ageList1}{
%%     ~~Describe \code{ageList1} here~~
}
  \item{ageList2}{
%%     ~~Describe \code{ageList2} here~~
}
  \item{ageList3}{
%%     ~~Describe \code{ageList3} here~~
}
  \item{ageList5}{
%%     ~~Describe \code{ageList5} here~~
}
  \item{preferedAgeModels1}{
%%     ~~Describe \code{preferedAgeModels1} here~~
}
  \item{preferedAgeModels2}{
%%     ~~Describe \code{preferedAgeModels2} here~~
}
  \item{preferedAgeModels3}{
%%     ~~Describe \code{preferedAgeModels3} here~~
}
  \item{preferedAgeModels4}{
%%     ~~Describe \code{preferedAgeModels4} here~~
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
function (ageList1, ageList2, ageList3, ageList5, preferedAgeModels1, 
    preferedAgeModels2, preferedAgeModels3, preferedAgeModels4) 
{
    blankPlot = ggplot(peakDF1) + geom_blank() + xlim(0, 10) + 
        ylim(0, 100)
    peakDF1 = data.frame(age = ageList1)
    ageHist1 <- ggplot(peakDF1) + geom_bar(binwidth = 1, aes(x = age, 
        fill = as.factor(1))) + geom_vline(xintercept = xmin) + 
        geom_vline(xintercept = xmax)
    peakDF2 = data.frame(age = ageList2[, 1], order = ageList2[, 
        2])
    ageHist2 <- ggplot(peakDF2) + geom_bar(binwidth = 1, aes(x = age, 
        fill = as.factor(order))) + geom_vline(xintercept = xmin) + 
        geom_vline(xintercept = xmax)
    peakDF3 = data.frame(age = ageList3[, 1], order = ageList3[, 
        2])
    ageHist3 <- ggplot(peakDF3) + geom_bar(binwidth = 1, aes(x = age, 
        fill = as.factor(order))) + geom_vline(xintercept = xmin) + 
        geom_vline(xintercept = xmax)
    peakDF4 = data.frame(age = ageList4[, 1], order = ageList4[, 
        2])
    ageHist4 <- ggplot(peakDF4) + geom_bar(binwidth = 1, aes(x = age, 
        fill = as.factor(order))) + geom_vline(xintercept = xmin) + 
        geom_vline(xintercept = xmax)
    if (length(preferedAgeModels1 > 0)) {
        peakDF1 = data.frame(age = preferedAgeModels1)
        chosenAgeHist1 <- ggplot(peakDF1) + geom_bar(binwidth = 1, 
            aes(x = age, fill = as.factor(1))) + geom_vline(xintercept = xmin) + 
            geom_vline(xintercept = xmax)
    }
    else {
        chosenAgeHist1 = blankPlot
    }
    if (length(preferedAgeModels2 > 0)) {
        peakDF2 = data.frame(age = preferedAgeModels2[, 1], order = preferedAgeModels2[, 
            2])
        chosenAgeHist2 <- ggplot(peakDF2) + geom_bar(binwidth = 1, 
            aes(x = age, fill = as.factor(order))) + geom_vline(xintercept = xmin) + 
            geom_vline(xintercept = xmax)
    }
    else {
        chosenAgeHist2 = blankPlot
    }
    if (length(preferedAgeModels3 > 0)) {
        peakDF3 = data.frame(age = preferedAgeModels3[, 1], order = preferedAgeModels3[, 
            2])
        chosenAgeHist3 <- ggplot(peakDF3) + geom_bar(binwidth = 1, 
            aes(x = age, fill = as.factor(order))) + geom_vline(xintercept = xmin) + 
            geom_vline(xintercept = xmax)
    }
    else {
        chosenAgeHist3 = blankPlot
    }
    if (length(preferedAgeModels4 > 0)) {
        peakDF4 = data.frame(age = preferedAgeModels4[, 1], order = preferedAgeModels4[, 
            2])
        chosenAgeHist4 <- ggplot(peakDF4) + geom_bar(binwidth = 1, 
            aes(x = age, fill = as.factor(order))) + geom_vline(xintercept = xmin) + 
            geom_vline(xintercept = xmax)
    }
    else {
        chosenAgeHist4 = blankPlot
    }
    theme_set(theme_bw(10))
    layout <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)
    multiplot(ageHist1, ageHist2, ageHist3, ageHist4, chosenAgeHist1, 
        chosenAgeHist2, chosenAgeHist3, chosenAgeHist4, layout = layout)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
