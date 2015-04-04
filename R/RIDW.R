#' Residual IDW interpolation
#'
#' @author Bastian Manz
#'
#' @description  This function calculates the residual (additive bias) at every 
#' pixel-point pair and interpolates them to the satellite field using inverse-
#' distance weighting (IDW) at every time step 
#'
#' @param gauge is a list: 
#' gauge[["points"]] is a data.frame with dimensions nrow=no of stations, 
#' containing columns x=x-coordinates, y=y-coordinates 
#' gauge[["ts"]] is a zoo object with dimensions ncol=no of stations, 
#' nrow=no of timestep
#' @param sat is a list
#' sat[["pixels"]] is a data.frame with dimensions nrow=no of satellite
#' pixels, containing columns x=x-coordinates, y=y-coordinates
#' sat  [["ts"]] is a zoo object with dimensions ncol=no of satellite 
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
#' # RIDW_out <- RIDW(sat, gauge)
#' # RIDW_cv  <- RIDW(sat, gauge, cross.val=TRUE)
#'

RIDW <- function(sat, gauge, longlat=TRUE, cross.val=FALSE){

Zs 	<- sat[[1]]
Zg 	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

# normalize vgm

#x-min(x)/min(x)-max(x)

# ensure Gaussian distribution (log-normal transformation, Box-Cox,normal score transform)

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata,Gdata[i,] 
                                                       ,longlat))

Zg       <- as.data.frame(t(Zg))
Zs_field <- as.data.frame(t(Zs))
Zs_trend <- as.data.frame(t(Zs[,loc]))


names(Zg)       <- gsub("*-","_",paste("Gauge_",names(Zs_trend),sep=""))
names(Zs_trend) <- gsub("*-","_",paste("Trend_",names(Zs_trend),sep=""))
names(Zs_field) <- names(Zs_trend)

gaugename<- names(Zg)
trendname<- names(Zs_trend)

## merge maps
data_df  <- cbind(Gdata, Zg, Zs_trend) 
data <- data_df

coordinates(data)     <- coordinates(Tdata[loc,])  
proj4string(data)     <- proj4string(Tdata)

coordinates(Zs_field) <- coordinates(Tdata) 
proj4string(Zs_field) <- proj4string(Tdata)

if (cross.val==FALSE){
  
  #resids_g  <- matrix(ncol=ncol(Zs_trend),nrow=nrow(Zs_trend))
  resids    <- matrix(ncol=ncol(Zs),nrow=nrow(Zs))
  Zs        <- matrix(ncol=nrow(Zs_field),nrow=length(gaugename))

  for (i in 1:length(gaugename)){
  
  print(i)
  
  # Get data for time step and exclude gauges with missing data
  data_sub1   <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]
  
  #resids_g[,i]<- Zg[,i] - Zs_trend[,i]
  resids_g <- Zg[,i] - Zs_trend[,i]
  resids_g <- resids_g[which(is.finite(resids_g))] 
  
  data_sub1@data[,gaugename[i]]  <- resids_g#[,i]
  formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
  resids[i,] <- idw(formula,locations=data_sub1,newdata=Zs_field,idp=2)$var1.pred 
  Zs[i,] <- as.numeric(as.matrix(Zs_field[,i]@data)) + resids[i,]
  
  Zs[i,Zs[i,]<0] =0

  }


  # tag dates 
       Zs  <- as.zoo(Zs); time(Zs) <- time(gauge[[1]])

  return(Zs) 

} else if(cross.val==TRUE){
  
  resids_g    <- matrix(ncol=ncol(Zs_trend),nrow=(nrow(Zs_trend)-1))
  resids       <- matrix(ncol=ncol(Zs),nrow=nrow(Zs))
  crossval        <-  matrix(ncol=nrow(Zs_trend),nrow=ncol(Zs_trend)) 
  
  for(i in 1:length(gaugename)){
    
    # Get data for time step and exclude gauges with missing data
    data_sub2    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]
    data_sub2 <- data_sub2[,c(gaugename[i],trendname[i])]
    
    for(j in 1:nrow(Zg)){
      
      print(paste(i,".",j,sep=""))
      
      data_sub2x <- data_sub2[-j,]
      
      resids_g[,i]<- data_sub2x@data[,1] - data_sub2x@data[,2]
      data_sub2x@data[,1]  <- resids_g[,i]
      formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
      resids[i,] <- idw(formula,locations=data_sub2x,newdata=Zs_field,idp=2)$var1.pred 
      crossval[i,j] <- as.numeric(as.matrix(Zs_field[,i]@data))[loc[j]] + resids[i,loc[j]]
 
    }
    
    crossval[i,crossval[i,]<0] =0
  }
  
  # tag dates 
  crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])
  
  return(crossval)
  
}
}

