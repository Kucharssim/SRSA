##########################################################
########### DEFINE USEFUL FUNCTIONS FOR SRSA #############
##########################################################

# Define the update function of the SR matrix
update <- function(M, I, a, g, i, j) { a * (I[,j] + g*M[,j] - M[,i]) }

# SRSA function - D=data for one participant on one trial
SRSA <- function(D, nAOI, a, g, mat=TRUE, M=matrix(0, nAOI, nAOI)){
  I <- diag(nAOI)
  for(t in 1:nrow(D)){
    i <- D[t,1]
    j <- D[t,2]
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
        Average.Matrices[which(i==Participants),] <- rowMeans(Matrices[which(i==Participants),,])
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
              model=lmod,
              projections=projections,
              predictions=predictions,
              strategies=pca$vectors[,1:nStrat],
              weightedstrategy=weighted.strategy)
  )
}