\name{crossval_score}
\alias{crossval_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
crossval_score(results_CV, gauge)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results_CV}{
%%     ~~Describe \code{results_CV} here~~
}
  \item{gauge}{
%%     ~~Describe \code{gauge} here~~
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
function (results_CV, gauge) 
{
    scorelist <- c("me", "mae", "r", "rmse", "nse")
    n <- length(results_CV)
    m <- ncol(gauge[[1]])
    o <- length(scorelist)
    subscore <- rep(NA, n * m * o)
    dim(subscore) <- c(m, n, o)
    for (i in 1:n) {
        for (j in 1:ncol(gauge[[1]])) {
            sim <- results_CV[[i]][, j]
            obs <- gauge[[1]][, j]
            subscore[j, i, 1] <- me(sim, obs, na.rm = TRUE)
            subscore[j, i, 2] <- mae(sim, obs, na.rm = TRUE)
            subscore[j, i, 4] <- rmse(sim, obs, na.rm = TRUE)
            subscore[j, i, 3] <- cor(sim, obs, "pairwise.complete.obs")
            subscore[j, i, 5] <- NSE(sim, obs, na.rm = TRUE)
        }
    }
    scores <- list()
    for (o in 1:length(scorelist)) {
        tmp <- data.frame(subscore[, , o])
        scores[[o]] <- tmp
    }
    names(scores) <- scorelist
    return(scores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
