\name{RIDW}
\alias{RIDW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RIDW(sat, gauge, longlat = TRUE, cross.val = FALSE)
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
    loc <- numeric()
    for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata, 
        Gdata[i, ], longlat))
    Zg <- as.data.frame(t(Zg))
    Zs_field <- as.data.frame(t(Zs))
    Zs_trend <- as.data.frame(t(Zs[, loc]))
    names(Zg) <- gsub("*-", "_", paste("Gauge_", names(Zs_trend), 
        sep = ""))
    names(Zs_trend) <- gsub("*-", "_", paste("Trend_", names(Zs_trend), 
        sep = ""))
    names(Zs_field) <- names(Zs_trend)
    gaugename <- names(Zg)
    trendname <- names(Zs_trend)
    data_df <- cbind(Gdata, Zg, Zs_trend)
    data <- data_df
    coordinates(data) <- coordinates(Tdata[loc, ])
    proj4string(data) <- proj4string(Tdata)
    coordinates(Zs_field) <- coordinates(Tdata)
    proj4string(Zs_field) <- proj4string(Tdata)
    if (cross.val == FALSE) {
        resids <- matrix(ncol = ncol(Zs), nrow = nrow(Zs))
        Zs <- matrix(ncol = nrow(Zs_field), nrow = length(gaugename))
        for (i in 1:length(gaugename)) {
            print(i)
            data_sub1 <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])), 
                ]
            resids_g <- Zg[, i] - Zs_trend[, i]
            resids_g <- resids_g[which(is.finite(resids_g))]
            data_sub1@data[, gaugename[i]] <- resids_g
            formula <- as.formula(paste(as.character(gaugename[i]), 
                " ~ ", 1, sep = ""))
            resids[i, ] <- idw(formula, locations = data_sub1, 
                newdata = Zs_field, idp = 2)$var1.pred
            Zs[i, ] <- as.numeric(as.matrix(Zs_field[, i]@data)) + 
                resids[i, ]
            Zs[i, Zs[i, ] < 0] = 0
        }
        Zs <- as.zoo(Zs)
        time(Zs) <- time(gauge[[1]])
        return(Zs)
    }
    else if (cross.val == TRUE) {
        resids_g <- matrix(ncol = ncol(Zs_trend), nrow = (nrow(Zs_trend) - 
            1))
        resids <- matrix(ncol = ncol(Zs), nrow = nrow(Zs))
        crossval <- matrix(ncol = nrow(Zs_trend), nrow = ncol(Zs_trend))
        for (i in 1:length(gaugename)) {
            data_sub2 <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])), 
                ]
            data_sub2 <- data_sub2[, c(gaugename[i], trendname[i])]
            for (j in 1:nrow(Zg)) {
                print(paste(i, ".", j, sep = ""))
                data_sub2x <- data_sub2[-j, ]
                resids_g[, i] <- data_sub2x@data[, 1] - data_sub2x@data[, 
                  2]
                data_sub2x@data[, 1] <- resids_g[, i]
                formula <- as.formula(paste(as.character(gaugename[i]), 
                  " ~ ", 1, sep = ""))
                resids[i, ] <- idw(formula, locations = data_sub2x, 
                  newdata = Zs_field, idp = 2)$var1.pred
                crossval[i, j] <- as.numeric(as.matrix(Zs_field[, 
                  i]@data))[loc[j]] + resids[i, loc[j]]
            }
            crossval[i, crossval[i, ] < 0] = 0
        }
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
