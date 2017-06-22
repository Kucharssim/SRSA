##############################################
################# SIMULATION #################
##############################################


########## Top-to-bottom Strategy ############

# Helping functions

assignProb1 <- function(probs, df){
  df <- 1-df
  probs <- probs / (sum(probs))
  probs * df
}

assignMat1 <- function(noise=0){
  require(MASS)
  require(arm)
  diagmeans <- c(0.62, 0.41, 0.45, 0.51, 0.27, 0.66)
  
  diagcov <- matrix(0.18, ncol=6, nrow=6) 
  diag(diagcov) <- 0.4 
  
  matTrans <- matrix(0, nrow=6, ncol=6)
  
  diag(matTrans) <- mvrnorm(1, diagmeans, diagcov)
  diag(matTrans) <- invlogit(diag(matTrans))
  
  matTrans[c(2,6),1] <- assignProb1(c(0.2, 0.1), matTrans[1,1])
  matTrans[c(1,3,6),2] <- assignProb1(c(0.2, 0.2, 0.1), matTrans[2,2])
  matTrans[c(2,4,6),3] <- assignProb1(c(0.1, 0.1, 0.2), matTrans[3,3])
  matTrans[c(3,5,6),4] <- assignProb1(c(0.1, 0.2, 0.1), matTrans[4,4])
  matTrans[c(1,3,4,6),5] <- assignProb1(c(0.1, 0.1, 0.14,0.1), matTrans[5,5])
  matTrans[1:3,6] <- assignProb1(c(0.15, 0.1, 0.1), matTrans[6,6])
  
  matTrans <- matTrans + runif(36, 0, noise)
  
  matTrans <- sweep(matTrans, 2, colSums(matTrans), "/")
  
  matTrans
}

fakeTrans1 <- function(){
  x <- rpois(1, 4.05)
  if(x %in% 0:1){
    x <- sample(3:4, 1, prob=c(2/3,1/3))
  }
  else if(x %in% c(5,6)){
    x <- sample(c(3:6,9:10), 1)
  }
  
  x
}

fakeEnd1 <- function(){
  if(rbinom(1,1,0.8)==1){
    x <- rpois(1, 3.1)
  }
  else { x <- 0}
  if(x==1){
    x <- sample(1:3,1, prob=c(0.35,0.35,0.3))
  }
  x
}

togglingMat1 <- function(noise){
  M <- matrix(c(0.65, 0.33, 0.05, 0.3 , 0.25,
                0.2 , 0.52, 0.1 , 0.05, 0.25,
                0.05, 0.05, 0.75, 0.2 , 0.125,
                0.05, 0.05, 0.05, 0.4 , 0.125,
                0.05, 0.05, 0.05, 0.05, 0.25), ncol=5, byrow=TRUE)
  M <- M + runif(25, 0, noise)
  M <- sweep(M, 2, colSums(M), "/")
  M
}


# The main function
simppt1 <- function(noise=0.05) {
  trans.to.six <- fakeTrans1()
  mat <- assignMat1(noise)
  toggling <- togglingMat1(noise)
  n.trans.to.six <- 0
  before <- 1
  
  path <- sample(1:6, 1, prob=c(0.57, 0.2, 0.01, 0.01, 0.01, 0.2))
  
  while (n.trans.to.six < trans.to.six | length(unique(path)) < 4) {
    sender <- path[length(path)]
    
    if (length(path) > 1) {
      if (sender == 6 & path[length(path) - 1] != 6) {
        before <- path[length(path) - 1]
      }
    }
    
    if (sender == 6 & length(path) != 1) {
      stay <- rbinom(1,1,mat[6,6])
      if(stay){
        receiver <- 6
      }
      else{
        receiver <- sample(1:5, 1, TRUE, toggling[, before])
      }
    }
    else{
      receiver <-  sample(1:6, 1, TRUE, mat[, sender])
    }
    
    path <- c(path, receiver)
    
    if (sender != 6 & receiver == 6) {
      n.trans.to.six <- n.trans.to.six + 1
    }
    
  }
  
  #add some sixes to the end
  c(path, rep(6, fakeEnd1()))
}






########## Systematic strategy #########

assignProb2 <- function(probs, df){
  df <- 1-df
  probs <- probs / (sum(probs))
  probs * df
}

assignMat2 <- function(noise=0){
  
  matTrans <- matrix(c(0.47, 0.03, 0.33, 0.12, 0, 0.06, 
                       0.15, 0.26, 0.36, 0.08, 0.03, 
                       0.12, 0.01, 0.02, 0.6, 0.02, 0, 
                       0.35, 0.03, 0.06, 0.47, 0.06, 0.25, 
                       0.12, 0, 0.29, 0.55, 0, 0.17, 
                       0, 0.01, 0.05, 0.51, 0.01, 0, 0.42), 
                     ncol=6,
                     byrow=FALSE)
  matTrans <- matTrans + runif(36, 0, noise)
  
  matTrans <- sweep(matTrans, 2, colSums(matTrans), "/")
  
  matTrans
}

fakeTrans2 <- function(){
  x <- rpois(1, 2.72)
  if(x == 0){
    x <- sample(2:3, 1, prob=c(2/3,1/3))
  }
  else if(x==1){
    x <- sample(1:3, 1, prob=c(1/5, 3/5, 1/5))
  }
  
  x
}

fakeEnd2 <- function(){
  if(rbinom(1,1,0.85)==1){
    x <- rpois(1, 2.3)
  }
  else { x <- 0}
  x
}

togglingMat2 <- function(noise){
  M <- matrix(c(0, 0  , 0   , 0, 1/5,
                0, 1/3, 0.05, 0, 1/5,
                1, 2/6, 0.95, 1, 1/5,
                0, 0  , 0   , 0, 1/5,
                0, 0  , 0   , 0, 1/5), ncol=5, byrow=TRUE)
  M <- M + runif(25, 0, noise)
  M <- sweep(M, 2, colSums(M), "/")
  M
}

simppt2 <- function(noise=0.05){
  trans.to.six <- fakeTrans2()
  mat <- assignMat2(noise)
  toggling <- togglingMat2(noise)
  n.trans.to.six <- 0
  before <- 1
  
  path <- sample(1:6, 1, prob=c(0.05, 0.20, 0.4, 0.05, 0.05, 0.3))
  
  
  while (n.trans.to.six < trans.to.six) {
    sender <- path[length(path)]
    
    if (length(path) > 1) {
      if (sender == 6 & path[length(path) - 1] != 6) {
        before <- path[length(path) - 1]
      }
    }
    
    if (sender == 6 & length(path) != 1) {
      stay <- rbinom(1,1,mat[6,6])
      if(stay){
        receiver <- 6
      }
      else{
        receiver <- sample(1:5, 1, TRUE, toggling[, before])
      }
    }
    else{
      receiver <-  sample(1:6, 1, TRUE, mat[, sender])
    }
    
    path <- c(path, receiver)
    
    if (sender != 6 & receiver == 6) {
      n.trans.to.six <- n.trans.to.six + 1
    }
    
  }
  
  #add some sixes to the end
  c(path, rep(6, fakeEnd2()))
}



########## Wrap under simulation #############

simulateFixed <- function(n, n.strategy.one, d, items=1, noise=0){
  strategy <- c(rep(1, n.strategy.one), rep(0, n-n.strategy.one))
  
  D <- data.frame(id=integer(0), item=integer(0), i=integer(0), j=integer(0))
  for (ppt in 1:n){
    for(i in 1:items){
      pattern <- vector()
      if(strategy[ppt]){pattern <- simppt1(noise)}
      else {pattern <- simppt2(noise)}
      
      D <- rbind(D, data.frame(id=ppt, item=i, i=pattern[-length(pattern)], j=pattern[-1]))
    }
  }
  
  Score <- ifelse(strategy, 
                  scale(rnorm(sum(strategy), 0, 1)),
                  scale(rnorm(n-sum(strategy), 0, 1))+d
                  )
  
  list(D=D, 
       ppt=data.frame(id=1:n, strategy=strategy, score=Score)
  )
}


repSimulateFixed <- function(k, n, n.strategy.one, d=1, items=1, noise=0){
  replicate(k, 
            simulateFixed(n, 
                          n.strategy.one, 
                          d, 
                          items=1, 
                          noise=0),
            simplify = FALSE)
}


