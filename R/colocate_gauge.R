#' Colocate gauges
#'
#' @author Zed Zulkafli
#'
#' @description This function colocates and average rain gauges estimates to the
#' nearest center cell of a target interpolation grid
#'
#' @param gauge is a list: 
#' gauge[[2]] is a data.frame with dimensions nrow=no of stations, 
#' containing columns x=x-coordinates, y=y-coordinates 
#' gauge[[1]] is a zoo object with dimensions ncol=no of stations, 
#' nrow=no of timestep
#' @param coords is sat  [[2]] is a data.frame with dimensions nrow=no of
#' satellite pixels, containing columns x=x-coordinates, y=y-coordinates
#' @param longlat is a flag for geospatial projection
#'
#' @return gauge see param description
#'
#' @examples
#' # gauge <- colocage.gauge(gauge, sat[[2]], longlat=TRUE)
#'

colocate.gauge <- function(gauge, coords, longlat=TRUE){

ts  <- gauge[[1]]
pts <- gauge[[2]]

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(pts)) loc[i] <- which.min(spDists(coords,pts[i,] ,longlat))

ts_colocated <- matrix(ncol = length(unique(loc)), nrow= nrow(ts))

for (t in 1: nrow(ts)) ts_colocated[t,] <- tapply(t(ts)[,t], loc, mean, na.rm=TRUE)

ts_colocated[!is.finite(ts_colocated)] <- NA

length_isna <- numeric()
for(i in 1:ncol(ts_colocated)) length_isna[i] <- 
  length(which(is.na(ts_colocated[,i])))

remove <- which(length_isna == nrow(ts_colocated))
if (length(remove > 1)) ts_colocated <- ts_colocated[,-remove]

ts_colocated <- as.zoo(ts_colocated); index(ts_colocated) <- index(ts); 
names(ts_colocated) <- 1:ncol(ts_colocated)

pts <- coords[unique(loc)[order(unique(loc))],]
if (length(remove > 1)) pts <- pts[-remove,]
pts$ID <- 1:nrow(pts)  
  
gauge[[1]] <- ts_colocated
gauge[[2]] <- pts

return(gauge) 

}