#' Bayesian Combination 
#'
#' @author Zed Zulkafli, Daniele Nerini
#'
#' @description This function merges satellite and Ordinary-Kriged rain 
#' gauge fields using Bayesian inference
#' References: Todini, E. (2001). A Bayesian technique for 
#' conditioning radar precipitation estimates to rain-gauge measurements. 
#' Hydrology and Earth System Sciences, 5(2):187-199. 
#' PROGEA Srl (2009). RAINMUSIC, User manual & references. PROGEA Srl, Bologna.
#'
#' #' @param gauge is a list: 
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
#' # BC_out <- BC(sat, gauge)
#' # BC_cv  <- BC(sat, gauge, cross.val=TRUE)
#'
################################################################################

# BC

################################################################################
BC <- function(sat,gauge,longlat=TRUE,cross.val=FALSE){

#library(rgdal)
#library(zoo)
#library(gstat)
#library(automap)

if(!cross.val){

Zs 	<- sat[[1]]
Zg 	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

# first create gridded rain gauge using OK 

Zg         <- data.frame(t(Zg))
gaugename  <- names(Zg)


## log transform rain gauge data (comment out if not needed)
#Zg <- log(Zg +0.01)

## merge maps
data <- cbind(Gdata, Zg)
coordinates(data) <- coordinates(Gdata) #~coords.x1+coords.x2 
proj4string(data) <- proj4string(Gdata)

vm.fit    <-  crossval <- maps <- list()

for (i in 1:length(gaugename)){

# Get data for time step and exclude gauges with missing data
  data_sub    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]

# Model semivariogram
  formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
  vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
                                  model = c("Sph", "Exp", "Gau"))

# Perform Kriging 
  maps[[i]] <-  krige0(formula, data=data_sub, newdata=sat[[2]],
                       model= vm.fit[[i]]$var_model, computeVar=TRUE,
                       fullCovariance = TRUE)
	
}

# back-log transform output
#	maps@data  			<- exp(maps@data) - 0.01

################################################################################

vm.fit2 <- maps2 <- list()

# calculate distances for generating covariance matrix
distances <- spDists(Tdata,Tdata,longlat)

for (i in 1:length(gaugename)){

#get Kriging output
  G 	     <- maps[[i]][[1]] 
  varG     <- maps[[i]][[2]]
  varG [varG < 0] <- 0

  S     <- Zs[i,]

# Calculate error
	errS <- S - G; errS  <- data.frame(t(errS)); names(errS) <- "errS"
	coordinates(errS) <- coordinates(Tdata)
	proj4string(errS) <- proj4string(Tdata)

# Model error semivariogram
      	formula     <- as.formula("errS ~ 1")
      	vm.fit2[[i]]<- autofitVariogram(formula, errS, model="Exp")

# Generate semivariance matrix between all pixels, then convert to covariance
      	Veps        <- variogramLine(object=vm.fit2[[i]]$var_model, dist_vector=distances)
      	Veps        <- vm.fit2[[i]]$var_model[2,"psill"] - Veps
      	Veps [Veps < 0] <- 0


# Run Kalman Filter

# A priori - Pstar
	mu_S  <- mean(t(S)[G>0] - G[G>0])
	Pstar <- t(S) - mu_S

# Kalman gain
	K  <- solve((varG+Veps),Veps) 

# Innovation
	Nu <- G - Pstar

# A posteriori
	BC <- Pstar + K %*% Nu       # prior at grid + correlation*observation
	BC[BC<0] <- 0
	
	Zs[i,] <- BC
}

        Zs  <- as.zoo(Zs); time(Zs) <- time(gauge[[1]])

return(Zs)

} else {

  
################################################################################

# BC - cross validation

################################################################################

crossval  <- matrix(ncol=nrow( gauge[[2]]), nrow=nrow(gauge[[1]]))

# initiate progress bar
pb <- txtProgressBar()
print("Bayesian combination - cross validation")

for (p in 1:length(gauge[[2]])){
setTxtProgressBar(pb, p/nrow(gauge[[2]]))

Zs 	<- sat[[1]]
Zg 	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

loc_p   <- which.min(spDists(Tdata,Gdata[p,] ,longlat))

Zg    	<- Zg[,-p]
Gdata 	<- Gdata[-p,]

# first create gridded rain gauge using OK 

Zg         <- data.frame(t(Zg))
gaugename  <- names(Zg)

## log transform rain gauge data (comment out if not needed)
#Zg <- log(Zg +0.01)

## merge maps
data <- cbind(Gdata, Zg)
coordinates(data) <- coordinates(Gdata) #~coords.x1+coords.x2 
proj4string(data) <- proj4string(Gdata)

vm.fit     <-  maps <- list()

for (i in 1:length(gaugename)){

# Get data for time step and exclude gauges with missing data
	data_sub    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]

# Model semivariogram
      	formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,
                                        sep=""))
      	vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
                                        model = c("Sph", "Exp", "Gau"))

# Perform Kriging 
      	maps[[i]] <-  krige0(formula, data=data_sub, newdata=sat[[2]],
      	                     model= vm.fit[[i]]$var_model, computeVar=TRUE,
      	                     fullCovariance = TRUE)
}


# back-log transform output
#	maps@data  			<- exp(maps@data) - 0.01


################################################################################

vm.fit2 <- maps2 <- list()

# calculate distances for generating covariance matrix
distances <- spDists(Tdata,Tdata,longlat)

for (i in 1:length(gaugename)){

#get Kriging output
  G 	     <- maps[[i]][[1]] 
  varG     <- maps[[i]][[2]]
  varG [varG < 0] <- 0

  S     <- Zs[i,]

# Calculate error
	errS <- S - G; errS  <- data.frame(t(errS)); names(errS) <- "errS"
	coordinates(errS) <- coordinates(Tdata)
	proj4string(errS) <- proj4string(Tdata)

# Model error semivariogram
       	formula     <- as.formula("errS ~ 1")
       	vm.fit2[[i]]<- autofitVariogram(formula, errS, model="Exp")

# Generate semivariance matrix between all pixels, then convert to covariance
      	Veps        <- variogramLine(object=vm.fit2[[i]]$var_model, 
                                     dist_vector=distances)
      	Veps        <- vm.fit2[[i]]$var_model[2,"psill"] - Veps
      	Veps [Veps < 0] <- 0

# Run Kalman Filter

# A priori - Pstar
	mu_S  <- mean(t(S)[G>0] - G[G>0])
	Pstar <- t(S) - mu_S

# Kalman gain
  K  <- solve((varG+Veps),Veps)  

# Innovation
	Nu <- G - Pstar

# A posteriori
	BC <- Pstar + K %*% Nu       # prior at grid + correlation*observation
	BC[BC<0] <- 0
	
	crossval[i,p] <- BC[loc_p]
}

}
close(pb)

      	crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])


return(crossval)
}

}



