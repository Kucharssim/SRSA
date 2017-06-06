##### Alternative SRSA ######

SR <- function(d, nAOI, g, mat=FALSE){
  i <- factor(d[,1], levels = 1:nAOI)
  j <- factor(d[,2], levels = 1:nAOI)
  trans <- table(j,i)
  trans <- prop.table(trans, 2)
  trans[is.nan(trans)] <- 0
  sr <- trans %*% solve(diag(nAOI)-g*trans, tol=1e-17)
  
  if(mat){
    return(sr)
  } else {
    return(c(sr))
  }
}

comp.SR <- function(fullD, nAOI, g){
  ppt.D <- split(fullD, fullD[,1])

  SRs <- sapply(ppt.D, function(D){
    item.D <- split(D, D[,2])

    srs <- sapply(item.D, function(d){
      SR(d[,3:4], nAOI, g)
    })

    return(rowMeans(srs))
  })

  return(t(SRs))
}

