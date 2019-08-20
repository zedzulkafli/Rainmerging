#' Ordinary Kriging
#'
#' @author Zed Zulkafli
#'
#' @description This function interpolate rain gauge data into the target grid 
#' using Ordinary Kriging at every time step
#'
#' @param gauge is a list: 
#' gauge[[2]] is a data.frame with dimensions nrow=no of stations, 
#' containing columns x=x-coordinates, y=y-coordinates 
#' gauge[[1]] is a zoo object with dimensions ncol=no of stations, 
#' nrow=no of timestep
#' @param coords is a data.frame with dimensions nrow=no of target grid pixels, 
#' containing columns x=x-coordinates, y=y-coordinates
#' @param cross.val option TRUE=in cross validation mode; default FALSE 
#'
#' @return Zs merging output - zoo object with dimensions == dimensions sat
#' [["ts""]
#' @return crossval cross validation- zoo object with dimensions == dimensions 
#' gauge[["ts""]
#'
#' @examples
#' # OK_out <- OK(sat, gauge)
#' # OK_cv  <- OK(sat, gauge, cross.val=TRUE)
#'

OK  <- function(coords, gauge, cross.val=FALSE){

Zg 	<- gauge[[1]]
Gdata 	<- gauge[[2]]

Zg         <- data.frame(t(Zg))
gaugename  <- names(Zg)

## log transform rain gauge data (comment out if not needed)
#Zg <- log(Zg +0.01)

## merge maps
data <- cbind(data.frame(Gdata), Zg)
coordinates(data) <- ~coords.x1+coords.x2 
proj4string(data) <- proj4string(Gdata)

vm.fit    <- list()

Zs        <- matrix(ncol=nrow(coords),nrow=length(gaugename))
crossval  <- matrix(ncol=nrow(Gdata), nrow=length(gaugename))


for (i in 1:length(gaugename)){

# Get data for time step and exclude gauges with missing data
	data_sub    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]

# Model semivariogram
      formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,
                                      sep=""))
      vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
                                      model = c("Sph", "Exp", "Gau"))

# Perform Kriging 
      Zs[i,] <- krige(formula, locations=data_sub, newdata=coords, 
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
