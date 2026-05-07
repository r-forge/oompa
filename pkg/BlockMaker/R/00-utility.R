truncscale <- function(x, dir=c("all", "row","col", "none"), ci=0.9, K=NULL) {
  x <- as.matrix(x)
  dir <- match.arg(dir)
  nc <- ncol(x)
  nr <- nrow(x)
  if(dir=="row") temp <- t(scale(t(x)))
  else if(dir=="col") temp <- scale(x)
  else if(dir=="all") temp <- scale(as.vector(x))
  else temp <- x
  if(!is.null(K)) {
    bound <- c(-K, K)
  } else {
    bound <- quantile(temp, c((1-ci)/2, (1+ci)/2), na.rm=TRUE)
  }
  temp[temp < bound[1]] <- bound[1]
  temp[temp > bound[2]] <- bound[2]
  ifelse(dir=="all",
         return(matrix(temp, nr,nc, dimnames=list(rownames(x), colnames(x)))),
         return(temp))
}
