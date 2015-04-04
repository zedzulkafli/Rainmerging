\name{OK}
\alias{OK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OK(coords, gauge, cross.val = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{gauge}{
%%     ~~Describe \code{gauge} here~~
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
function (coords, gauge, cross.val = FALSE) 
{
    Zg <- gauge[[1]]
    Gdata <- gauge[[2]]
    Zg <- data.frame(t(Zg))
    gaugename <- names(Zg)
    data <- cbind(Gdata, Zg)
    coordinates(data) <- ~coords.x1 + coords.x2
    proj4string(data) <- proj4string(Gdata)
    vm.fit <- list()
    Zs <- matrix(ncol = nrow(coords), nrow = length(gaugename))
    crossval <- matrix(ncol = nrow(Gdata), nrow = length(gaugename))
    for (i in 1:length(gaugename)) {
        data_sub <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])), 
            ]
        formula <- as.formula(paste(as.character(gaugename[i]), 
            " ~ ", 1, sep = ""))
        vm.fit[[i]] <- autofitVariogram(formula, data_sub, model = c("Sph", 
            "Exp", "Gau"))
        Zs[i, ] <- krige(formula, locations = data_sub, newdata = coords, 
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
