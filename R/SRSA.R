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
