plotSR <- function(values, title=NULL, includeVal = FALSE,
                   nAOI=sqrt(length(values)), g=NULL){
  m <- matrix(values, ncol=nAOI)
  if(g==1) {stop("gamma cannot be equal to 1")}
  if(!is.null(g)) {
    minmax <- c(0, 1/(1-g))
  } else {
    minmax <- c(0, max(m))
  }
  corrplot(m, method="shade", is.corr = FALSE,
           tl.col = "black", tl.srt = 0, tl.offset = 0.7,
           p.mat=m, insig = ifelse(includeVal, "p-value", "n"),
           sig.level = min(m)-1, cl.lim=minmax)
  text(nAOI/2+0.5, nAOI + 1.5, title, cex=2)
}
