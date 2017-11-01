### New functions ###

srsa <- function(data, nAOI, method = 'pca', nStrat = 2, nComp = NULL, target,
                 resolution = 1/50, collapse = FALSE, converged = FALSE,
                 a = NULL, g = NULL){
  
  if(collapse){
    data <- data[data[,3] != data[,4],]
  }
  
  # iterating and not converged - supply starting values
  if(!is.null(a) & !is.null(g) & converged == FALSE){
    cat('Iterating... Searching optimal alpha and gamma \n')
    
    output <- DoSRSA(data, nAOI, target, a, g)
    class(output) <- 'srsaPCA'
    return(output)
  } else if(!is.null(g) & converged == TRUE){
    cat('Iterating... Searching optimal gamma \n')
    cat('Currently not working')
    
    output <- DoSRSA(data, nAOI, target, a, g)
    class(output) <- 'srsaPCA'
    return(output)
  }
  
  
  
  grid <- seq(0,1,resolution)
  grid <- expand.grid(a=grid, g=grid)
  
  if(method == 'pca'){
    
    
  } else if(method == 'kmeans'){
    
  } else {
    stop('The method is not recognized, current options are pca and kmeans')
  }
}

summary.srsaPCA <- function(object, ...){
  cat(object$par)
}

print.srsaPCA <- function(x, ...){
  cat(x$par)
}