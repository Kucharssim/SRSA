# Functions for SR matrix construction #

#' Create SR given one scanpath
#' @description This function takes two vectors of states (i and j) and converts it into Successor Representation of the transition sequence.
 
#' @param data Matrix or dataframe with two columns (first column is the sender (i), second receiver (j))
#' @param nAOI positive integer indicating the number of AOIs
#' @param a alpha. Numeric between 0 and 1
#' @param g gamma. Numeric between 0 and 1
#' @param converged Logical. Should the SR matrix be constructed statically based on first order transition matrix?
#' @param mat Logical. Should the output be matrix or vector? Used for internal purposes.
 
#' @return matrix or vector of the Successor Representation values. The vector is columnwise unfold of a matrix.

sr <- function(data, nAOI, a=0, g=0, converged = FALSE, mat=TRUE){
  data <- data[,3:4]
  I <- diag(nAOI) 
  M <- matrix(0, nAOI, nAOI)
  
  if(converged){
    if(g >= 1) {stop('gamma has to be less than 1.')}
    
    i <- factor(data$i, levels = 1:nAOI)
    j <- factor(data$j, levels = 1:nAOI)
    trans <- table(j,i)
    trans <- prop.table(trans, 2)
    M <- trans %*% solve(I - g * trans, tol = 1e-17)
    
  } else {
    for(t in seq_along(data[[1]])){
      i <- data$i[t]
      j <- data$j[t]
      M[,i] <- M[,i] + updateSR(M, I, a, g, i, j) 
    }
  }
  
  if(mat){
    return(M)
  } else {
    return(c(M))
  }
  
}

# Update function
updateSR <- function(M, I, a, g, i, j) { a * (I[,j] + g*M[,j] - M[,i]) }


###### Create SR matrices per participant ######

#' Create a Successor Representations
#' @param data matrix or dataframe of 4 columns. 1 Id of participant, 2 item, 3 sender and 4 receiver
#' @param nAOI number of Areas of interest (positive integer)
#' @param a alpha. Numeric between 0 and 1
#' @param g gamma. Numeric between 0 and 1
#' @param converged Logical. Should the SR matrix be constructed statically based on first order transition matrix?
#' @param averageItems Logical. Should the matrices be averaged across items?
#' @param normalize Logical. Should the matrix cells be normalized?
srs <- function(data, nAOI, a, g, converged = FALSE,
                averageItems = TRUE, normalize = FALSE, parallel=FALSE){

  # compute the SR matrix per scanpath
  M <-  plyr::ddply(.data=data, .variables = c('id', 'item'),
                    .fun = sr, 
                    nAOI=nAOI, a=a, g=g, mat=FALSE, 
                    .drop = FALSE, .parallel = parallel)
  
  if(averageItems) {
    
    M <- aggregate(M, list(M$id), mean, na.rm = TRUE)[-c(1:3)]

    if(normalize) {
      suppressWarnings(
        M <- scale(M)
      )
      M[is.nan(M)] <- 0
    }
  }
  
  return(as.matrix(M))
}
