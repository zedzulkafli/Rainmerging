#' Cross validation scores
#'
#' @author Zed Zulkafli
#'
#' @description This function calculates performance indices for 
#' for the cross-validation output for comparison between the methods
#'
#' @param gauge is a list: 
#' gauge[[2]] is a data.frame with dimensions nrow=no of stations, 
#' containing columns x=x-coordinates, y=y-coordinates 
#' gauge[[1]] is a zoo object with dimensions ncol=no of stations, 
#' nrow=no of timestep
#' @param results_CV is a list of the output from one or multiple merging functions run in cross-validation mode
#'
#' @return scores array containing performance indices for the cross-validation 
#'
#' @examples
#' # results <- crossval_score(results_CV,gauge)
#'

crossval_score <- function(results_CV,gauge){

#library(hydroGOF)

scorelist <- c("me", "mae",  "r", "rmse", "nse")

n <- length(results_CV)
m <- ncol(gauge[[1]])
o <- length(scorelist)

subscore <- rep(NA,n*m*o)
dim(subscore) <- c(m,n,o)


for (i in 1:n){

for (j in 1:ncol(gauge[[1]])){
sim 		    <- results_CV[[i]][,j]
obs 		    <- gauge[[1]][,j]
subscore[j,i,1]     <- me(sim,obs,na.rm=TRUE)
subscore[j,i,2]     <- mae(sim,obs,na.rm=TRUE)
subscore[j,i,4]     <- rmse(sim,obs,na.rm=TRUE)
if (length(sim)> 1) subscore[j,i,3]     <- cor(sim,obs,"pairwise.complete.obs")
if (length(sim)> 1) subscore[j,i,5]     <- NSE(sim,obs,na.rm=TRUE)
}
}


scores <- list()
for(o in 1:length(scorelist)) {
tmp 	    <- data.frame(subscore [,,o])
names(tmp)  <- names(results_CV)
scores [[o]]    <- tmp
}

names(scores) <- scorelist

return(scores)
}
