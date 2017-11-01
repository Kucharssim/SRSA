### New functions ###

srsa <- function(data, nAOI, method = 'pca', nStrat = 2, nComp = NULL, 
                 target = NULL, resolution = 1/50, collapse = FALSE,
                 converged = FALSE, sort.cor = TRUE, rotate = FALSE,
                 a = NULL, g = NULL, parallel = FALSE){
  
  if(collapse){
    data <- data[data[,3] != data[,4],]
  }
  
  # PCA analysis
  if(method=='pca'){
    
    if(!converged & is.null(g) & is.null(a)){
      cat('Approximating the parameter space... Resolution:', resolution)
      
      map <- mapPCA(data, nAOI, target, nStrat, nComp, sort.cor, parallel, rotate)
      output <- optimPCA(data, nAOI, target, map$par[1], map$par[2],
                         nStrat, nComp, sort.cor, rotate, parallel)
      
    } else if(!is.null(a) & !is.null(g) & !converged){
      cat('Iterating... Searching optimal alpha and gamma \n')
      
      output <- optimPCA(data, nAOI, target, a, g,
                         nStrat, nComp, sort.cor, rotate, parallel)
      
    } else if(!is.null(g) & converged){
      cat('Iterating... Searching optimal gamma \n')
      cat('Currently not working')
      
      output <- DoSRSA(data, nAOI, target, a, g)
      
    } else {
      stop('The function cannot be executed')
    }
    
    class(output) <- 'srsaPCA'
    return(output)
    
    
  } else if(method == 'kmeans'){ # k means analysis
    
    if(!converged & is.null(g) & is.null(a)){
      map <- mapKmeans(data, nAOI, target, grid, nStrat, parallel)
    
      output <- optimKmeans(data, nAOI, map$par[1], map$par[2], nStrat, parallel)

    } else if(!converged & !is.null(g) & !is.null(a)){
      
      output <- optimKmeans(data, nAOI, map$par[1], map$par[2], nStrat, parallel)
      
    } else  if(converged & !is.null(g)){
      # make the srs
      
      # compute the k-means
    } else {
      stop('The function cannot be executed')
    }
    
    class(output) <- 'srsaKmeans'
    return(output)
  }
}


#### cross validation ####
# 
# cv.srsaPCA(D, nAOI, target, nStrat=2, nComp=NULL, sort.cor=TRUE,
#            rotate=FALSE, parallel=FALSE, covariate=NULL, resolution){
#   
#   pb <- txtProgressBar(min = 0, max = length(criterion), style = 3)
#   fits <- lapply(unique(D[,1]), function(leave.out){
#     
#     train <- as.data.frame(D[D[,1]!=leave.out,])
#     test <- as.data.frame(D[D[,1]==leave.out,])
#     
#     leave <- which(unique(D[,1])==leave.out)
#     setTxtProgressBar(pb, leave)
#     
#     par <- mapPCA(train, nAOI, criterion[-leave], nStrat, nComp, covariate,
#                   resolution, rotate)
#     
#     res <- srsaPCA(train, nAOI, criterion[-leave],
#                   par$par[1], par$par[2], nStrat, nComp,
#                   sort.cor, rotate, parallel=FALSE,
#                   covariate)
#     
#     if(abs(par$par[1]-res$par[1])>resolution |
#        abs(par$par[2]-res$par[2])>resolution){
#       warning("the solution did not converge towards the global optimum. Participant:", leave.out)
#     }
#     
#     sr.leave <- comp.SRSA(test, 6, res$par[1], res$par[2], TRUE, FALSE)
#     mean.sr <- colMeans(res$rSRSA, na.rm = TRUE)
#     sd.sr <- apply(res$rSRSA, 2, sd, na.rm = TRUE)
#     
#     sr.leave <- (sr.leave-mean.sr)/sd.sr
#     sr.leave[abs(sr.leave)==Inf] <- 0
#     sr.leave[is.na(sr.leave)] <- 0
#     
#     pred <- sr.leave %*% res$weightedstrategy + res$coefficients[1,1]
#     list(prediction=pred, parameters=res$par)
#   })
#   
#   fits
#   
# }