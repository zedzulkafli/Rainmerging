#' Mean Bias Correction
#'
#' @author Zed Zulkafli
#'
#' @description  This function calculates a spatial mean multiplicative bias 
#' correction of the satellite field at every time step 
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
#' # MBC_out <- MBC(sat, gauge)
#' # MBC_cv  <- MBC(sat, gauge, cross.val=TRUE)
#'

MBC <- function(sat,gauge, longlat=TRUE, cross.val=FALSE){

Zs 	<- sat[[1]]
Zg 	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

if (!cross.val) {

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata,Gdata[i,] ,
                                                       longlat))

# subset Zs 
Zs_sub <- Zs[,loc]

# remove pixel-point pairs with NA values
Zs_sub[is.na(Zg)] <- NA

# calculate bias factor
CF   <- apply(Zg, 1, sum, na.rm=TRUE)/ apply(Zs_sub, 1, sum,  na.rm=TRUE) 

# apply bias correction 
Zs   <- CF*Zs

return( Zs ) 

} else {

################################################################################

# MBC  - cross validation

################################################################################

crossval <- matrix(ncol=nrow(Gdata),nrow=nrow(Zg))

# initiate progress bar
pb <- txtProgressBar()
print("Mean bias correction - cross validation")

for (p in 1:nrow(Gdata)){
setTxtProgressBar(pb, p/nrow(Gdata))

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(Gdata[-p])) loc[i] <- which.min(spDists(Tdata,Gdata[i,] ,
                                                           longlat=TRUE))

# subset Zs 
Zs_sub <- Zs[,loc]

# remove pixel-point pairs with NA values
Zs_sub[is.na(Zg[,-p])] <- NA

# calculate bias factor
CF   <- apply(Zg[,-p], 1, sum, na.rm=TRUE)/ apply(Zs_sub, 1, sum,  na.rm=TRUE) 

# apply bias correction 
loc_p          <- which.min(spDists(Tdata,Gdata[p,] ,longlat=TRUE))
crossval[,p]   <- CF*Zs[,loc_p]
}
close(pb)

crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])

return( crossval)
}

}
