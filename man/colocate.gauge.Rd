\name{colocate.gauge}
\alias{colocate.gauge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
colocate.gauge(gauge, coords, longlat = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gauge}{
%%     ~~Describe \code{gauge} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{longlat}{
%%     ~~Describe \code{longlat} here~~
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
function (gauge, coords, longlat = TRUE) 
{
    ts <- gauge[[1]]
    pts <- gauge[[2]]
    loc <- numeric()
    for (i in 1:length(pts)) loc[i] <- which.min(spDists(coords, 
        pts[i, ], longlat))
    ts_colocated <- matrix(ncol = length(unique(loc)), nrow = nrow(ts))
    for (t in 1:nrow(ts)) ts_colocated[t, ] <- tapply(t(ts)[, 
        t], loc, mean, na.rm = TRUE)
    ts_colocated[!is.finite(ts_colocated)] <- NA
    length_isna <- numeric()
    for (i in 1:ncol(ts_colocated)) length_isna[i] <- length(which(is.na(ts_colocated[, 
        i])))
    remove <- which(length_isna == nrow(ts_colocated))
    if (length(remove > 1)) 
        ts_colocated <- ts_colocated[, -remove]
    ts_colocated <- as.zoo(ts_colocated)
    index(ts_colocated) <- index(ts)
    names(ts_colocated) <- 1:ncol(ts_colocated)
    pts <- coords[unique(loc)[order(unique(loc))], ]
    if (length(remove > 1)) 
        pts <- pts[-remove, ]
    pts$ID <- 1:nrow(pts)
    gauge[[1]] <- ts_colocated
    gauge[[2]] <- pts
    return(gauge)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
