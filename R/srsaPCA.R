#### main PCA function ####
srsaPCA <- function(data, nAOI, target, a, g, nStrat, nComp,
                    sort.cor=TRUE, converged = FALSE, rotate=FALSE,
                    parallel=FALSE, covariate = NULL, optim='both'){
  
  par <- optimPCA(data, nAOI, target, a, g, nStrat, nComp,
                  sort.cor, converged, rotate, parallel=FALSE, 
                  covariate, optim)
  
  rmatrices <- srs(data, nAOI, par[1], par[2], FALSE, TRUE, FALSE)
  matrices <- suppressWarnings(scale(rmatrices))
  matrices[is.na(matrices)] <- 0
  #matrices <- srs(data, nAOI, par[1], par[2], FALSE, TRUE, TRUE)
  
  pca <- projectPCA(matrices, target, nStrat, nComp, sort.cor, rotate)
  
  if(is.null(covariate)){
    model <- lm(target ~ pca$projections)
  } else {
    predictors <- cbind(pca$projections, covariate)
    model <- lm(target ~ predictors)
  }
  
  strategies <- pca$pca$vectors[, pca$comp.order[1:nStrat]]
  weightMat <- as.matrix(strategies) %*% t(t(model$coefficients[2:(nStrat+1)]))
  
  return(list(par=par,
              matrices=matrices,
              rmatrices=rmatrices,
              PCA=pca,
              coefficients=summary(model)$coefficients,
              r.squared=summary(model)$r.squared,
              #projections=projections,
              predictions=predict(model),
              strategies=strategies,
              weights=weightMat)
  )
}


optimPCA <- function(data, nAOI, target, a, g, nStrat, nComp=NULL,
                     sort.cor = TRUE, converged = FALSE,
                     rotate = FALSE, parallel = FALSE,
                     covariate = NULL, optim = 'both'){
  if(optim=='both'){
    o <- optim(c(a, g), costPCA, gr=NULL, 
               data, nAOI, target, nStrat, nComp,
               sort.cor, converged, rotate,
               covariate, optim, NULL, 
               method="L-BFGS-B", lower=0, upper=1)$par
    
  } else if(optim=='a'){
    a <- optim(a, costPCA, gr=NULL,
               data, nAOI, target, nStrat, nComp,
               sort.cor, converged, rotate,
               covariate, optim, g,
               method="L-BFGS-B", lower=0, upper=1)$par
    o <- c(a, g)
  } else {
    g <- optim(g, costPCA, gr=NULL,
               data, nAOI, target, nStrat, nComp,
               sort.cor, converged, rotate,
               covariate, optim, a,
               method="L-BFGS-B", lower=0, upper=1)$par
    o <- c(a, g)
  }
  
  return(o)
}

#### Predict ####
costPCA <- function(par, data, nAOI, target, nStrat, nComp,
                    sort.cor=TRUE, converged = FALSE, 
                    rotate=FALSE, covariate=NULL,
                    optim='both', placeholder){
  
  if(optim=='both'){
    matrices <- srs(data, nAOI, par[1], par[2], converged, TRUE, TRUE)
  } else if(optim=='a') {
    matrices <- srs(data, nAOI, par, placeholder, converged, TRUE, TRUE)
  } else if(optim=='g') {
    matrices <- srs(data, nAOI, placeholder, par, converged, TRUE, TRUE)
  } else {
    stop('Something went wrong: Cannot optimize with the current settings')
  }
  
  pca <- projectPCA(matrices, target, nStrat, nComp, sort.cor, rotate)
  
  if(is.null(covariate)){
    model <- lm(target ~ pca$projections)
  } else {
    predictors <- cbind(pca$projections, covariate)
    model <- lm(target ~ predictors)
  }
  
  return(1 - summary(model)$r.squared)
}

projectPCA <- function(matrices, target, nStrat, nComp,
                       sort.cor=TRUE, rotate=FALSE){
  suppressWarnings(
    corM <- cor(matrices)
  )
  corM[is.na(corM)] <- 0
  
  pca <- eigen(corM)
  
  # if(rotate=='varimax'){
  #   
  # } else if (rotate=='promax'){
  #   
  # }
  
  if(is.null(nComp)){
    nComp <- length(pca$values)
  }
  
  if(sort.cor){
    projections <- matrices %*% pca$vectors
    projections <- projections[, 1:nComp, drop=FALSE]
    comp.order <- orderCor(projections, target)
    
    projections <- projections[,comp.order][,1:nStrat, drop=FALSE]
  } else {
    projections <- matrices %*% pca$vectors[,1:nStrat, drop=FALSE]
    comp.order <- 1:nComp
  }
  
  return(list(projections=projections,
              pca = pca,
              comp.order = comp.order)
  )
}

orderCor <- function(projections, target){
  suppressWarnings(
    correlations <- abs(apply(projections, 2, cor, y = target))
  )
  
  return(order(-correlations))
}

# map PCA R squared

mapPCA <- function(data, nAOI, target, resolution, nStrat, nComp=NULL,
                   sort.cor = TRUE, converged = FALSE,
                   rotate = FALSE, parallel = FALSE,
                   covariate = NULL, optim = 'both', placeholder=0){
  
  grid <- seq(0, 1, resolution)
  
  if(optim=='both'){
    values <- grid
    grid <- expand.grid(a=grid, g=grid)
    map <- mapply(function(a, g) {costPCA(c(a, g), data, nAOI, target, 
                                          nStrat, nComp, sort.cor, converged,
                                          rotate, covariate, optim)}, 
                  grid$a, grid$g, SIMPLIFY = FALSE)
    map <- t(matrix(1 - unlist(map), ncol=sqrt(nrow(grid))))
    
    colnames(map) <- rownames(map) <- values
    par <- values[which(map==max(map), arr.ind=TRUE)]
    
  } else if(optim %in% c('a', 'g')){
    map <- sapply(grid, costPCA, data, nAOI, target, nStrat, nComp, sort.cor,
                  converged, rotate, covariate, optim, placeholder)
    map <- 1 - map
    names(map) <- grid
    par <- grid[which.max(map)]
    
    if(optim=='a') {
      par <- c(par, placeholder)
    } else {
      par <- c(placeholder, par)  
    }
    
  } else {
    stop('The settings is not valid (mapPCA)')
  }
  
  return(list(map=map,
              par=par)
         )
}


### Cross-validation ###

# cv.srsaPCA(data, nAOI, target, resolution, nStrat, nComp=NULL,
#            sort.cor = TRUE, converged = FALSE,
#            rotate = FALSE, parallel = FALSE,
#            covariate = NULL, optim = 'both', placeholder=0,
#            detailed = FALSE){
#   
#   id <- unique(data[, 1])
#   
#   folds <- lapply(id, crossvalidate.srsaPCA,
#                   data=data, nAOI=nAOI, target=target, resolution=resolution,
#                   nStrat=nStrat, nComp=nComp, sort.cor=sort.cor, 
#                   converged=converged, rotate=rotate, parallel=parallel,
#                   covariate=covariate, optim=optim, placeholder=placeholder,
#                   detailed=detailed)
#   
# }
# 
# crossvalidate.srsaPCA <- function(leave, data, nAOI, target, resolution,
#                                   nStrat, nComp=NULL,
#                                   sort.cor = TRUE, converged = FALSE,
#                                   rotate = FALSE, parallel = FALSE,
#                                   covariate = NULL, optim = 'both',
#                                   placeholder=0, detailed = FALSE){
#   
#   which.leave <- which(D[,1]!=leave)
#   train <- as.data.frame(D[D[,1]!=leave,])
#   train.target <- target[-which.leave]
#   
#   test <- as.data.frame(D[D[,1]==leave,])
#   test.targe <- target[which.leave]
#   
#   map <- mapPCA(train, nAOI, train.target, resolution, nStrat, nComp, 
#                 sort.cor, converged, rotate, parallel, covariate, optim, 
#                 placeholder)
#   
# }