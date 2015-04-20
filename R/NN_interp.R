#' Nearest-neighbour interpolator
#'
#' @author Zed Zulkafli
#'
#' @description This function interpolates satellite data to a new grid using 
#' the nearest-neighbour interpolation
#'
#' @param sat is a list
#' sat[[2]] is a data.frame with dimensions nrow=no of satellite
#' pixels, containing columns x=x-coordinates, y=y-coordinates, ID=pixel ID
#' sat  [[1]] is a zoo object with dimensions ncol=no of satellite 
#' pixels, nrow=no of timestep
#' @param coords is a data.frame with dimensions nrow=no of target grid pixels, 
#' containing columns x=x-coordinates, y=y-coordinates
#'
#' @return sat see param description
#'
#' @examples
#' # sat <- NN.interp(sat, coords)
#'


NN.interp <- function(sat, coords){

#read in source data coordinates
map     <- sat[[2]]
ts 	    <- sat[[1]]

#interpolate source data ID to interpolation grid 
map 	<- idw(ID ~ 1, map, coords, nmax = 1)[,"var1.pred"]

#interpolate entire time series
ts      <- as.data.frame(ts)[,as.character(map$var1.pred)]
ts      <- as.zoo(ts); index(ts) <- index(sat[[1]])

#pack output
sat[[1]]   <- ts
sat[[2]]   <- map
  
return(sat)
}