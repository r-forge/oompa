# Copyright (C) Kevin R. Coombes, 2007-2016

##########################################################
# Graphical utilities

fliplr <-
# MATLAB clone to flip matrices
function(x) {
  return(x[, seq(ncol(x), 1, -1)])
}

flipud <-
# MATLAB clone to flip matrices
function(x)
{
  return(x[seq(nrow(x), 1, -1),  ])
}

ellipse <-
# Adds an ellipse to a plot, with axes oriented along the coordinate axes.
function(a, b, x0=0, y0=0, ...) {
  theta <- (0:360)*3.14159/180
  x <- x0[1] + a*cos(theta)
  y <- y0[1] + b*sin(theta)
  points(x, y, ...)
}

f.qq <- 
# just a normal QQ-plot, with some extra lines.
function(x, main='', cut=0, ...) {
  qqnorm(x, ylab='Quantiles of Test Statistic', main=main, ...)
  qqline(x)
  if (cut) {
    abline(h=c(cut, -cut))
    text(-2, cut+0.6, paste('cutoff:', cut), adj=0)
  }
  invisible(x)
}

f.qt <-
# a QQ-plot against the t-distribution, with some extra lines.
function(x, df, main='', cut=0, ...) {
  y <- qt(ppoints(length(x)), df)
  data.quartiles <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  norm.quartiles <- qt(c(0.25, 0.75), df)
  b <- (data.quartiles[2] - data.quartiles[1])/(norm.quartiles[2] -
                                                norm.quartiles[1])
  a <- data.quartiles[1] - norm.quartiles[1] * b
  
  qqplot(y, x, ylab='Quantiles of Test Statistic',
         xlab=paste('Quantiles of T distribution, df =', df), main=main, ...)
  abline(a, b, ...)
  if (cut) {
    abline(h=c(cut, -cut))
    text(-2, cut+0.6, paste('cutoff:', cut), adj=0)
  }
  invisible(x)
}

##########################################################################
# Miscellaneous utilities

f.above.thresh <-
# Returns the fraction of elements in the vector A that are greater than 
# the threshold t.
function(a, t = 1)
{
  return(length(a[a > t])/length(a))
}

f.cord <-
# Computes and returns the concordance coefficient
# of the two input vectors x and y.
function(x, y, inf.rm=FALSE)
{
# straightforward computation of concordance
  if(inf.rm) {
    x[abs(x)==Inf] <- NA
    y[abs(y)==Inf] <- NA
  }
  bad.points <- is.na(x) | is.na(y)
  if (sum(bad.points) > 0) {	# need this for cross-platform compatability
                  # R and S-Plus use different arguments to var to omit NA's.
    x <- x[!bad.points]
    y <- y[!bad.points]
  }
  vxy <- var(x, y)
  vx <- var(x)
  vy <- var(y)
  mx <- mean(x,na.rm=TRUE)
  my <- mean(y,na.rm=TRUE)
  return((2 * vxy)/(vx + vy + (mx - my)^2))
}

f.oneway.rankings <- function(r, s) { order(s)[r] }
