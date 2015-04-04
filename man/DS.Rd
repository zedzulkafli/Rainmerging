\name{DS}
\alias{DS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DS(sat, gauge, longlat = TRUE, cross.val = FALSE)
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
    index(Zs) <- as.Date(index(Zs))
    index(Zg) <- as.Date(index(Zg))
    if (!cross.val) {
        loc <- numeric()
        for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata, 
            Gdata[i, ], longlat))
        Zs_sub <- Zs[, loc]
        Zs_sub[is.na(Zg)] <- NA
        res <- Zs_sub - Zg
        distres <- spDists(Gdata, Tdata, longlat)
        distpseud <- spDists(Tdata, Tdata, longlat)
        kerf <- function(x, b) {
            return(1/sqrt(2 * pi) * exp(-x^2/(2 * b^2)))
        }
        b1 <- 1.06 * sd(distres) * length(distres)^(-1/5)
        b2 <- 1.06 * sd(distpseud) * length(distpseud)^(-1/5)
        K1 <- kerf((distres), b1)
        K2 <- kerf((distpseud), b2)
        K1 <- t(K1)
        K1 <- rep(K1, nrow(res))
        dim(K1) <- c(dim(t(distres)), nrow(res))
        for (t in 1:nrow(res)) {
            out <- which(is.na(res[t, ]))
            K1[, out, t] <- 0
        }
        res[which(is.na(res))] <- 0
        eDS <- eSS <- matrix(ncol = nrow(Tdata), nrow = nrow(res))
        for (t in 1:nrow(res)) {
            eSS[t, ] <- (K1[, , t] \%*\% t(res[t, ]))/apply(K1[, 
                , t], 1, sum)
            eDS[t, ] <- (K1[, , t] \%*\% t(res[t, ]) + K2 \%*\% eSS[t, 
                ])/(apply(K1[, , t], 1, sum) + apply(K2, 1, sum))
        }
        Zs <- Zs - eDS
        Zs[which(Zs < 0)] <- 0
        return(Zs)
    }
    else {
        crossval <- matrix(ncol = nrow(Gdata), nrow = nrow(Zg))
        pb <- txtProgressBar()
        print("Double kernel smoothing - cross validation")
        for (q in 1:nrow(Gdata)) {
            setTxtProgressBar(pb, q/nrow(Gdata))
            loc <- numeric()
            for (i in 1:length(Gdata[-q, ])) loc[i] <- which.min(spDists(Tdata, 
                Gdata[i, ], longlat))
            Zs_sub <- Zs[, loc]
            Zs_sub[is.na(Zg[, -q])] <- NA
            res <- Zs_sub - Zg[, -q]
            distres <- spDists(Gdata[-q, ], Tdata, longlat)
            distpseud <- spDists(Tdata, Tdata, longlat)
            kerf <- function(x, b) {
                return(1/sqrt(2 * pi) * exp(-x^2/(2 * b^2)))
            }
            b1 <- 1.06 * sd(distres) * length(distres)^(-1/5)
            b2 <- 1.06 * sd(distpseud) * length(distpseud)^(-1/5)
            K1 <- kerf((distres), b1)
            K2 <- kerf((distpseud), b2)
            K1 <- t(K1)
            K1 <- rep(K1, nrow(res))
            dim(K1) <- c(dim(t(distres)), nrow(res))
            for (t in 1:nrow(res)) {
                out <- which(is.na(res[t, ]))
                K1[, out, t] <- 0
            }
            res[which(is.na(res))] <- 0
            eSS <- matrix(ncol = nrow(Tdata), nrow = nrow(res))
            eDS <- numeric()
            loc_q <- which.min(spDists(Tdata, Gdata[q, ], longlat))
            for (t in 1:nrow(res)) {
                eSS[t, ] <- (K1[, , t] \%*\% t(res[t, ]))/apply(K1[, 
                  , t], 1, sum)
                eDS[t] <- (t(matrix(K1[loc_q, , t])) \%*\% t(res[t, 
                  ]) + K2[, loc_q] \%*\% eSS[t, ])/(sum(K1[loc_q, 
                  , t]) + sum(K2[, loc_q]))
            }
            DS <- Zs[, loc_q] - eDS
            DS[which(DS < 0)] <- 0
            crossval[, q] <- DS
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
