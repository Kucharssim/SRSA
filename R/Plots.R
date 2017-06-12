plotSR <- function(values, title=NULL, includeVal = FALSE,
                   nAOI=sqrt(length(values))){
  m <- matrix(values, ncol=nAOI)
  corrplot(m, method="shade", is.corr = FALSE,
           tl.col = "black", tl.srt = 0, tl.offset = 0.7,
           p.mat=m, insig = ifelse(includeVal, "p-value", "n"),
           sig.level = min(m)-1)
  text(nAOI/2+0.5, nAOI + 1.5, title, cex=2)
}
