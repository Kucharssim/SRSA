# code for k-means srsa

srsaKmeans <- function(){}

optimKmeans <- function(){}

mapKmeans <- function(data, nAOI, resolution = 1/100, nStrat, replicates = 1,
                      converged = FALSE, parallel = FALSE, normalize = FALSE,
                      optim = 'both', fixed.par = 0){
  
  grid <- seq(resolution, 1-resolution, by = resolution)
  if(nStrat == 1){
    return(1)
  }
  
  if(optim=="both"){
    values <- grid
    grid <- expand.grid(a=grid, g=grid)
    
    map <- mapply(function(a, g){
      screeKmeans(a, g, data, nAOI, nStrat, replicates, converged, normalize)
    }, grid$a, grid$g, SIMPLIFY = FALSE)
    
    map <- do.call(rbind, map)
    colnames(map) <- 2:nStrat
    map <- cbind(grid, '1' = 1, map)
    
  } else if(optim=="a"){
    map <- sapply(grid, screeKmeans, 
                  fixed.par, data, nAOI, nStrat, replicates, converged, normalize)
    
    #colnames(map)<- 2:nStrat
    map <- cbind(a = grid, g = fixed.par, '1' = 1, map)
    
  } else if(optim=="g"){
    map <- sapply(grid, function(g){
      screeKmeans(fixed.par, g, data, nAOI, nStrat, replicates, converged, normalize)
    })
    
    #colnames(map)<- 2:nStrat
    map <- cbind(a = fixed.par, g = grid, '1' = 1, map)
    
    }
  
  return(map)
}

costKmeans <- function(k, matrices, replicates){
  k.means <- replicate(replicates, kmeans(matrices, k), simplify = FALSE)
  
  totwith <- sapply(k.means, function(x){
    x$tot.withinss / x$totss
  })
  
  uniques <- unique(totwith)
  if(length(uniques) > 1) {
    warning('There are multiple local optima (k-means): ',
            k, ' clusters, values: ', 
            paste(round(uniques, 2)), collapse = ", ")
  }
  
  return(list(
    value = min(totwith),
    uniques = uniques
  ))
}

screeKmeans <- function(a, g, data, nAOI, nStrat, replicates = 1,
                        converged = FALSE, normalize = FALSE){
  if(nStrat != 1){
    matrices <- srs(data, nAOI, a, g, converged, normalize = normalize)
    res <- lapply(2:nStrat, costKmeans, matrices, replicates)
  
    totwith <- sapply(res, function(r) r$value)
  } else {
    totwith <- NA
  }
  
  return(totwith)
}

