\name{KED}
\alias{KED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KED(sat, gauge, longlat = TRUE, cross.val = FALSE)
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
    Zg <- data.frame(t(Zg))
    Zs_field <- data.frame(t(Zs))
    Zs_trend <- data.frame(t(Zs[, loc]))
    names(Zg) <- paste("Gauge", names(Zg), sep = "")
    names(Zs_trend) <- paste("Trend", names(Zs_trend), sep = "")
    names(Zs_field) <- names(Zs_trend)
    gaugename <- names(Zg)
    trendname <- names(Zs_trend)
    data <- cbind(Gdata, Zg, Zs_trend)
    coordinates(data) <- coordinates(Tdata[loc, ])
    proj4string(data) <- proj4string(Tdata)
    coordinates(Zs_field) <- coordinates(Tdata)
    proj4string(Zs_field) <- proj4string(Tdata)
    vm.fit <- list()
    crossval <- matrix(ncol = nrow(Gdata), nrow = length(gaugename))
    for (i in 1:length(gaugename)) {
        data_sub <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])), 
            ]
        formula <- as.formula(paste(as.character(gaugename[i]), 
            " ~ ", as.character(trendname[i]), sep = ""))
        vm.fit[[i]] <- autofitVariogram(formula, data_sub, model = c("Sph", 
            "Exp", "Gau"))
        Zs[i, ] <- krige(formula, locations = data_sub, newdata = Zs_field, 
            model = vm.fit[[i]]$var_model)$var1.pred
        if (cross.val == TRUE) {
            cv <- krige.cv(formula, locations = data_sub, model = vm.fit[[i]]$var_model)
            crossval[i, as.numeric(row.names(cv))] <- cv$var1.pred
        }
    }
    Zs <- as.zoo(Zs)
    time(Zs) <- time(gauge[[1]])
    crossval <- as.zoo(crossval)
    time(crossval) <- time(gauge[[1]])
    if (cross.val == TRUE) 
        return(crossval)
    else return(Zs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
