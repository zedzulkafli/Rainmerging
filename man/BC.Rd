\name{BC}
\alias{BC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BC(sat, gauge, longlat = TRUE, cross.val = FALSE)
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
    if (!cross.val) {
        Zs <- sat[[1]]
        Zg <- gauge[[1]]
        Tdata <- sat[[2]]
        Gdata <- gauge[[2]]
        Zg <- data.frame(t(Zg))
        gaugename <- names(Zg)
        data <- cbind(Gdata, Zg)
        coordinates(data) <- coordinates(Gdata)
        proj4string(data) <- proj4string(Gdata)
        vm.fit <- crossval <- maps <- list()
        for (i in 1:length(gaugename)) {
            data_sub <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])), 
                ]
            formula <- as.formula(paste(as.character(gaugename[i]), 
                " ~ ", 1, sep = ""))
            vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
                model = c("Sph", "Exp", "Gau"))
            maps[[i]] <- krige(formula, locations = data_sub, 
                newdata = sat[["pixels"]], model = vm.fit[[i]]$var_model)
        }
        vm.fit2 <- maps2 <- list()
        distances <- spDists(Tdata, Tdata, longlat)
        for (i in 1:length(gaugename)) {
            G <- maps[[i]]$var1.pred
            varG <- variogramLine(object = vm.fit[[i]]$var_model, 
                dist_vector = distances)
            varG <- vm.fit[[i]]$var_model[2, "psill"] - varG
            varG[varG < 0] <- 0
            S <- Zs[i, ]
            errS <- S - G
            errS <- data.frame(t(errS))
            names(errS) <- "errS"
            coordinates(errS) <- coordinates(Tdata)
            proj4string(errS) <- proj4string(Tdata)
            formula <- as.formula("errS ~ 1")
            vm.fit2[[i]] <- autofitVariogram(formula, errS, model = "Exp")
            Veps <- variogramLine(object = vm.fit2[[i]]$var_model, 
                dist_vector = distances)
            Veps <- vm.fit2[[i]]$var_model[2, "psill"] - Veps
            Veps[Veps < 0] <- 0
            mu_S <- mean(t(S)[G > 0] - G[G > 0])
            Pstar <- t(S) - mu_S
            K <- solve((varG + Veps), Veps)
            Nu <- G - Pstar
            BC <- Pstar + K \%*\% Nu
            BC[BC < 0] <- 0
            Zs[i, ] <- BC
        }
        Zs <- as.zoo(Zs)
        time(Zs) <- time(gauge[[1]])
        return(Zs)
    }
    else {
        crossval <- matrix(ncol = nrow(gauge[[2]]), nrow = nrow(gauge[[1]]))
        pb <- txtProgressBar()
        print("Bayesian combination - cross validation")
        for (p in 1:length(gauge[[2]])) {
            setTxtProgressBar(pb, p/nrow(gauge[[2]]))
            Zs <- sat[[1]]
            Zg <- gauge[[1]]
            Tdata <- sat[[2]]
            Gdata <- gauge[[2]]
            loc_p <- which.min(spDists(Tdata, Gdata[p, ], longlat))
            Zg <- Zg[, -p]
            Gdata <- Gdata[-p, ]
            Zg <- data.frame(t(Zg))
            gaugename <- names(Zg)
            data <- cbind(Gdata, Zg)
            coordinates(data) <- coordinates(Gdata)
            proj4string(data) <- proj4string(Gdata)
            vm.fit <- maps <- list()
            for (i in 1:length(gaugename)) {
                data_sub <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])), 
                  ]
                formula <- as.formula(paste(as.character(gaugename[i]), 
                  " ~ ", 1, sep = ""))
                vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
                  model = c("Sph", "Exp", "Gau"))
                maps[[i]] <- krige(formula, locations = data_sub, 
                  newdata = sat[["pixels"]], model = vm.fit[[i]]$var_model, 
                  debug.level = FALSE)
            }
            vm.fit2 <- maps2 <- list()
            distances <- spDists(Tdata, Tdata, longlat)
            for (i in 1:length(gaugename)) {
                G <- maps[[i]]$var1.pred
                varG <- variogramLine(object = vm.fit[[i]]$var_model, 
                  dist_vector = distances)
                varG <- vm.fit[[i]]$var_model[2, "psill"] - varG
                varG[varG < 0] <- 0
                S <- Zs[i, ]
                errS <- S - G
                errS <- data.frame(t(errS))
                names(errS) <- "errS"
                coordinates(errS) <- coordinates(Tdata)
                proj4string(errS) <- proj4string(Tdata)
                formula <- as.formula("errS ~ 1")
                vm.fit2[[i]] <- autofitVariogram(formula, errS, 
                  model = "Exp")
                Veps <- variogramLine(object = vm.fit2[[i]]$var_model, 
                  dist_vector = distances)
                Veps <- vm.fit2[[i]]$var_model[2, "psill"] - 
                  Veps
                Veps[Veps < 0] <- 0
                mu_S <- mean(t(S)[G > 0] - G[G > 0])
                Pstar <- t(S) - mu_S
                K <- solve((varG + Veps), Veps)
                Nu <- G - Pstar
                BC <- Pstar + K \%*\% Nu
                BC[BC < 0] <- 0
                crossval[i, p] <- BC[loc_p]
            }
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
