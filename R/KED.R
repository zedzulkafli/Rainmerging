#' Kriging with External Drift
#'
#' @author Zed Zulkafli, Bastian Manz
#'
#' @description This function performs a Kriging interpolation of rain gauge 
#' estimates using satellite estimates as an external drift
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
#' @return crossval cross validation- zoo object with dimensions == dimensions gauge
#' [["ts""]
#'
#' @examples
#' # KED_out <- KED(sat, gauge)
#' # KED_cv  <- KED(sat, gauge, cross.val=TRUE)
#'

KED <- function(sat,gauge, longlat =TRUE, cross.val=FALSE){

Zs 	<- sat[[1]]
Zg 	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

# get location of Zg in Zs

loc <- numeric()
for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata,Gdata[i,],
                                                       longlat))

Zg       <- data.frame(t(Zg))
Zs_field <- data.frame(t(Zs))
Zs_trend <- data.frame(t(Zs[,loc]))

names(Zg)       <- paste("Gauge",names(Zg),sep="")
names(Zs_trend) <- paste("Trend",names(Zs_trend),sep="")
names(Zs_field) <- names(Zs_trend)

gaugename<- names(Zg)
trendname<- names(Zs_trend)
 
## log transform rain gauge data (comment out if not needed)
#Zg <- log(Zg +0.01)

## merge maps
data  <- cbind(Gdata, Zg, Zs_trend)

coordinates(data) <- coordinates(Tdata[loc,]) 
proj4string(data)     <- proj4string(Tdata)

coordinates(Zs_field) <- coordinates(Tdata) 
proj4string(Zs_field) <- proj4string(Tdata)

vm.fit <-  list()
crossval  <- matrix(ncol=nrow(Gdata),nrow=length(gaugename))

for (i in 1:length(gaugename)){

# Get data for time step and exclude gauges with missing data
	data_sub    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]

# Model semivariogram
  formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 
                                  as.character(trendname[i]),sep=""))
  vm.fit[[i]] <- autofitVariogram(formula, data_sub,  
                                  model = c("Sph", "Exp", "Gau"))

# Perform Kriging
      Zs[i,] <- krige(formula, locations=data_sub, newdata=Zs_field, 
                      model = vm.fit[[i]]$var_model)$var1.pred 
	
# Perform cross-validation 	
      if(cross.val==TRUE) {
      cv <- krige.cv(formula, locations=data_sub, model = vm.fit[[i]]$var_model)
      crossval[i,as.numeric(row.names(cv))] <- cv$var1.pred
      }

}

# back-log transform output
# 	Zs    <- exp(Zs) - 0.01

# tag dates 
       Zs  <- as.zoo(Zs); time(Zs) <- time(gauge[[1]])
       crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])

# return results	   
if(cross.val==TRUE) return(crossval) else return(Zs) 
}

