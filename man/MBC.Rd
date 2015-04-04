\name{MBC}
\alias{MBC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MBC(sat, gauge, longlat = TRUE, cross.val = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sat}{
%%     ~~Describe \code{sat} here~~
}
  \item{gauge}{
%%     ~~Describe \code{gauge} here~~
}
  \item{longlat}{
%%     ~~Describe \code{longlat} here~~
}
  \item{cross.val}{
%%     ~~Describe \code{cross.val} here~~
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
function (sat, gauge, longlat = TRUE, cross.val = FALSE) 
{
    Zs <- sat[[1]]
    Zg <- gauge[[1]]
    Tdata <- sat[[2]]
    Gdata <- gauge[[2]]
    if (!cross.val) {
        loc <- numeric()
        for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata, 
            Gdata[i, ], longlat))
        Zs_sub <- Zs[, loc]
        Zs_sub[is.na(Zg)] <- NA
        CF <- apply(Zg, 1, sum, na.rm = TRUE)/apply(Zs_sub, 1, 
            sum, na.rm = TRUE)
        Zs <- CF * Zs
        return(Zs)
    }
    else {
        crossval <- matrix(ncol = nrow(Gdata), nrow = nrow(Zg))
        pb <- txtProgressBar()
        print("Mean bias correction - cross validation")
        for (p in 1:nrow(Gdata)) {
            setTxtProgressBar(pb, p/nrow(Gdata))
            loc <- numeric()
            for (i in 1:length(Gdata[-p])) loc[i] <- which.min(spDists(Tdata, 
                Gdata[i, ], longlat = TRUE))
            Zs_sub <- Zs[, loc]
            Zs_sub[is.na(Zg[, -p])] <- NA
            CF <- apply(Zg[, -p], 1, sum, na.rm = TRUE)/apply(Zs_sub, 
                1, sum, na.rm = TRUE)
            loc_p <- which.min(spDists(Tdata, Gdata[p, ], longlat = TRUE))
            crossval[, p] <- CF * Zs[, loc_p]
        }
        close(pb)
        crossval <- as.zoo(crossval)
        time(crossval) <- time(gauge[[1]])
        return(crossval)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
