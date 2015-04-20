#' Double kernel residual Smoothing
#'
#' @author Zed Zulkafli, Daniele Nerini
#'
#' @description This function merges satellite and rain gauge fields using 
#' the double kernel smoothing of the residuals following Li & #' Shao (2010)
#' References: Li, M. and Shao, Q. (2010). 
#' An improved statistical approach to merge satellite rainfall
#' estimates and raingauge data. Journal of Hydrology, 385(1-4):51-64.
#' Silverman, B. (1986). Density Estimation for Statistics and Data Analysis. 
#' Chapman and Hall/CRC Monographs on Statistics and Applied Probability Series. 
#' Chapman and Hall/CRC.
#'
#' @param gauge is a list: 
#' gauge[[2]] is a data.frame with dimensions nrow=no of stations, 
#' containing columns x=x-coordinates, y=y-coordinates 
#' gauge[[1]] is a zoo object with dimensions ncol=no of stations, 
#' nrow=no of timestep
#' @param sat is a list
#' sat[[2]] is a data.frame with dimensions nrow=no of satellite
#' pixels, containing columns x=x-coordinates, y=y-coordinates
#' sat  [[1]] is a zoo object with dimensions ncol=no of satellite 
#' pixels, nrow=no of timestep
#' @param cross.val option TRUE=in cross validation mode; default FALSE 
#' @param longlat is a flag to describe the coordinate grids of spatial data. 
#' Defaults to TRUE i.e. on long-lat notation 
#'
#' @return Zs merging output - zoo object with dimensions == dimensions sat 
#' [["ts""]
#' @return crossval cross validation- zoo object with dimensions == dimensions 
#' gauge[["ts""]
#'
#' @examples
#' # DS_out <- DS(sat, gauge)
#' # DS_cv  <- DS(sat, gauge, cross.val=TRUE)
#'
################################################################################

DS <- function(sat,gauge,longlat=TRUE, cross.val=FALSE){

Zs	<- sat[[1]]
Zg	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

# standardise index classes
index(Zs) <- as.Date(index(Zs))
index(Zg) <- as.Date(index(Zg))

if(!cross.val) {

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata,Gdata[i,] 
                                                       ,longlat))

# subset Zs 
Zs_sub <- Zs[,loc]

# get point residuals
res <- Zs_sub - Zg

# calculate kernel weights

# first level -  get distances between interpolation grids and point residuals
distres <- spDists(Gdata, Tdata ,longlat)

# second level - get distances between interpolation grids and pseudo-residuals
distpseud <- spDists(Tdata, Tdata ,longlat)

# Define Gaussian kernel function
kerf <- function(x,b){ return(1/sqrt(2*pi)*exp(- x^2 /(2*b^2))) } 

# Parameterise b using Silverman's rule of thumb

b1 <- 1.06*sd(distres  )*length(distres  )^(-1/5)
b2 <- 1.06*sd(distpseud)*length(distpseud)^(-1/5);

# compute Kernel weights 
K1 <- kerf((distres)  ,b1);
K2 <- kerf((distpseud),b2);

K1 <- t(K1)
K1 <- rep(K1,nrow(res)); dim(K1) <- c(dim(t(distres)),nrow(res))

# interpolate residuals
eDS <- eSS <- matrix(ncol=nrow(Tdata), nrow=nrow(res))

for (t in 1:nrow(res)){

eSS[t,] <- (K1[,,t] %*% t(res[t,])) /apply(K1[,,t],1,sum)
eDS[t,] <- (K1[,,t] %*% t(res[t,]) + K2 %*% eSS[t,] )/ (apply(K1[,,t],1,sum)
                                                        + apply(K2,1,sum))
}

# add residual back to background
Zs <- Zs - eDS

Zs[which(Zs<0)] <- 0 

return(Zs)
} else {

################################################################################

# DS  - cross validation

################################################################################

crossval <- matrix(ncol=nrow(Gdata),nrow=nrow(Zg))

# initiate progress bar
pb <- txtProgressBar()
print("Double kernel smoothing - cross validation")

for (q in 1:nrow(Gdata)){
setTxtProgressBar(pb, q/nrow(Gdata))

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(Gdata[-q,])) loc[i] <- which.min(spDists(Tdata,Gdata[i,] 
                                                            ,longlat))

# subset Zs 
Zs_sub <- Zs[,loc]

# get point residuals
res <- Zs_sub - Zg[,-q]

# calculate kernel weights

# first level -  get distances between interpolation grids and point residuals
distres   <- spDists(Gdata[-q,], Tdata ,longlat)

# second level - get distances between interpolation grids and pseudo-residuals
distpseud <- spDists(Tdata, Tdata ,longlat)

# Define Gaussian kernel function
kerf <- function(x,b){ return(1/sqrt(2*pi)*exp(- x^2 /(2*b^2))) } 

# Parameterise b using Silverman's rule of thumb
b1 <- 1.06*sd(distres  )*length(distres  )^(-1/5)
b2 <- 1.06*sd(distpseud)*length(distpseud)^(-1/5);
    
# compute Kernel weights 
K1 <- kerf((distres)  ,b1);
K2 <- kerf((distpseud),b2);

K1 <- t(K1)
K1 <- rep(K1,nrow(res)); dim(K1) <- c(dim(t(distres)),nrow(res))

# Interpolate residuals
eSS <- matrix(ncol=nrow(Tdata), nrow=nrow(res))
eDS <- numeric()

loc_q          <- which.min(spDists(Tdata,Gdata[q,] ,longlat))

for (t in 1:nrow(res)){

eSS[t,] <- (K1[,,t] %*% t(res[t,])) /apply(K1[,,t],1,sum)
eDS[t]  <- (t(matrix(K1[loc_q,,t])) %*% t(res[t,]) + K2[,loc_q] %*% eSS[t,] )/ 
             (sum(K1[loc_q,,t]) + sum(K2[,loc_q]))

}

DS 		 <- Zs[,loc_q] - eDS
DS[which(DS<0)] <- 0 

crossval[,q]    <- DS

}

close(pb)
crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])

return(crossval)
}
}
