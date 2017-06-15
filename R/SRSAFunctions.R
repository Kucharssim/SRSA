##########################################################
########### DEFINE USEFUL FUNCTIONS FOR SRSA #############
##########################################################
#
# Define the update function of the SR matrix
update <- function(M, I, a, g, i, j) { a * (I[,j] + g*M[,j] - M[,i]) }

# SRSA function - D=data for one participant on one trial
SRSA <- function(D, nAOI, a, g, mat=TRUE, M=matrix(0, nAOI, nAOI)){
  I <- diag(nAOI)
  
  if(nrow(D)==0){
    return(rep(NA, nAOI^2))
  }
  for(t in 1:nrow(D)){
    i <- D[t,1]
    j <- D[t,2]
    #print(1/t)
    M[,i] <- M[,i] + update(M, I, a, g, i, j)
  }
  
  if(mat){M}
  else(c(M))
}

# Compute SRSA for full dataset - col 1=ID, col 2=Item, col 3 and 4=sender and reciever
comp.SRSA <- function(fullD, nAOI, a, g, averageItems=TRUE, normalize=TRUE) {
  Participants <- unique(fullD[,1])
  Items <- unique(fullD[,2])
  
  #3 dim array: participants x matrix in a vector format x items
  Matrices <- replicate (length(Items), matrix(0, length(Participants), nAOI^2))
  
  for(k in Items){
    # subset of "k" item
    D_I <- fullD[fullD[,2]==k,]
    
    for(l in Participants){
      # subset of participant and just i and j columns
      D <- D_I[D_I[,1]==l,3:4]
      
      # Compute SR "matrix" (as a vector) and store it for the specific participant
      # on the specific item
      Matrices[which(l==Participants),,which(k==Items)] <- SRSA(D, nAOI, a, g, mat=FALSE)
    }
  }
  
  
  if(!averageItems) {
    # normalize the features
    if(normalize){
      Matrices <- apply(Matrices, c(2,3), scale)
    }
    if(length(Items)==1){
      Matrices <- Matrices[,,1]
    }
    
    return(Matrices)
  }
  
  # Average across items
  else{ 
    # participants x SR matrix
    Average.Matrices <- matrix(NA, dim(Matrices)[1], dim(Matrices)[2])
    
    if(length(Items)==1) {
      warning("Cannot average across 1 Item")
      Average.Matrices <- Matrices[,,1]}
    
    else{
      for (i in Participants){
        Average.Matrices[which(i==Participants),] <- rowMeans(Matrices[which(i==Participants),,], na.rm = TRUE)
      }
    }

    if(normalize){
      Average.Matrices <- apply(Average.Matrices, 2, scale)
    }
    return(Average.Matrices)
  }
}

# Cost Function:
# Input: SR Matrix obtained with current a and g, criterion what we want to predict,
# m components we want to use for prediction, k estimated components(default to what Hayes 2011 used)
# Estimates the components and predicts the score for each participant
# Correlates the score with criterion and select the two with the strongest one
# Use them in linear model predicting criterion
# OUTPUT: (1 - R Squared)
Cost <- function(pcaMatrix, criterion, nComp=2, hayes){
  pcor <- cor(pcaMatrix, use="p")
  pcor[is.na(pcor)] <- 0
  pca <- eigen(pcor)
  pcaMatrix[is.nan(pcaMatrix)] <- 0
  pca$vectors <- highestCorr(pca$vectors, criterion, pcaMatrix, hayes)
  projections <- pcaMatrix %*% pca$vectors[,1:nComp]
  
  lmod <- lm (criterion~projections)
  1 - summary(lmod)$r.squared
}

# Wrapper for Cost function (so it can be called with optim() )
# Currently only for one Item option
CostSRSA <- function(parameters, fullD, nAOI, criterion, nStrat, hayes){
  pcaMatrix <- comp.SRSA(fullD, nAOI, parameters[1], parameters[2])
  
  Cost(pcaMatrix, criterion, nStrat, hayes)
}

highestCorr <- function(vectors, criterion, data, hayes){
  if(hayes){
    projections <- data %*% vectors
    corrs <- abs(apply(projections, 2, cor, y=criterion))
    return(vectors[,order(-corrs)])
  }
  else {
    return(vectors)
  }
}

DoSRSA <- function(D, nAOI, criterion, alpha, gamma, nStrat=2, hayes=TRUE){
  o <- optim(c(alpha, gamma), CostSRSA, gr=NULL, D, nAOI, criterion, nStrat, hayes, 
             method="Nelder-Mead", lower=0, upper=1)
  
  matSRSA <- comp.SRSA(D, nAOI, o$par[1], o$par[2])
  unstdSRSA <- comp.SRSA(D, nAOI, o$par[1], o$par[2], normalize = FALSE)
  #matSRSA <- apply(unstdSRSA,2,scale)
  
  pcor <- cor(matSRSA, use="p")
  pcor[is.na(pcor)] <- 0
  pca <- eigen(pcor)
  matSRSA[is.nan(matSRSA)] <- 0
  pca$vectors <- highestCorr(pca$vectors, criterion, matSRSA, hayes)
  projections <- matSRSA %*% pca$vectors[,1:nStrat]
  
  lmod <- lm (criterion~projections)
  predictions <- lmod$coefficients[1] + projections %*% lmod$coefficients[2:(nStrat+1)]
  
  weighted.strategy <- pca$vectors[,1:nStrat] %*% t(t(lmod$coefficients[2:(nStrat+1)]))
  
  return(list(par=o$par,
              SRSA=matSRSA,
              rSRSA=unstdSRSA,
              PCA=pca,
              coefficients=summary(lmod)$coefficients,
              r.squared=summary(lmod)$r.squared,
              projections=projections,
              predictions=predictions,
              strategies=pca$vectors[,1:nStrat],
              weightedstrategy=weighted.strategy)
  )
}



optimR <- function(D, criterion, nStrat, hayes, resolution=1/100){
  D <- as.data.frame(D)
  criterion <- as.vector(criterion)
  require(parallel)
  require(SRSA)
  
  cores <- detectCores()-1
  cl <- makeCluster(cores)
  grid <- grid <- seq(0,1,by=resolution)
  e.grid <- expand.grid(a=grid, g=grid)
  clusterExport(cl, varlist=c("CostSRSA", "Cost", "SRSA",
                              "update", "comp.SRSA", "highestCorr"),
              envir = globalenv())
  R <- clusterMap(cl, function(a, g) {CostSRSA(c(a,g), D, 6, criterion, 2, TRUE)},
                  e.grid$a, e.grid$g)

  stopCluster(cl)
  R <- t(matrix(1 - unlist(R), ncol=length(grid)))
  colnames(R) <- rownames(R) <- grid

  filled.contour(grid, grid, R, zlim=c(0,ceiling(max(R)*10)/10), 
                 color.palette = colorRamps::matlab.like,
                 xlab="Gamma", ylab="Alpha", main="R^2 map")
  
  return(list(R=R[which(R==max(R), arr.ind = TRUE)],
              par=grid[which(R==max(R), arr.ind = TRUE)][2:1]))
}


crossSR <- function(D, nAOI, criterion, nStrat=2, hayes=TRUE){
  pb <- txtProgressBar(min = 0, max = length(criterion), style = 3)
  fits <- lapply(unique(D[,1]), function(leave.out){
    
    train <- as.data.frame(D[D[,1]!=leave.out,])
    test <- as.data.frame(D[D[,1]==leave.out,])
    
    leave <- which(unique(D[,1])==leave.out)
    setTxtProgressBar(pb, leave)
    
    #par(mfrow=rep(ceiling(sqrt(length(criterion))),2))
    par <- optimR(D, criterion, nStrat, hayes, 1/50)
    
    res <- DoSRSA(train, 6, criterion[-leave],
                  par$par[1], par$par[2], nStrat, hayes)
    
    if(abs(par$par[1]-res$par[1])>1/50 |
       abs(par$par[2]-res$par[2])>1/50){
      warning("the solution did not converge towards the global optimum. Participant:", leave.out)
    }
    
    sr.leave <- comp.SRSA(test, 6, res$par[1], res$par[2], TRUE, FALSE)
    mean.sr <- colMeans(res$rSRSA, na.rm = TRUE)
    sd.sr <- apply(res$rSRSA, 2, sd, na.rm = TRUE)
    
    sr.leave <- (sr.leave-mean.sr)/sd.sr
    sr.leave[abs(sr.leave)==Inf] <- 0
    sr.leave[is.na(sr.leave)] <- 0
    
    pred <- sr.leave %*% res$weightedstrategy
    list(prediction=pred, parameters=res$par)
  })
  
  fits
}
