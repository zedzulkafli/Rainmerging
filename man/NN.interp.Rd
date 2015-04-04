\name{NN.interp}
\alias{NN.interp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NN.interp(sat, coords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sat}{
%%     ~~Describe \code{sat} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
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
function (sat, coords) 
{
    map <- sat[[2]]
    ts <- sat[[1]]
    map <- idw(ID ~ 1, map, coords, nmax = 1)[, "var1.pred"]
    ts <- as.data.frame(ts)[, as.character(map$var1.pred)]
    ts <- as.zoo(ts)
    index(ts) <- index(sat[[1]])
    sat[[1]] <- ts
    sat[[2]] <- map
    return(sat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
