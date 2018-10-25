##############################
# internal computations

krc <- function(nj, yj, alpha=1, beta=1) {
  J <- length(nj)
  if (J != length(yj)) stop("bad sizes")
  temp <- - lgamma(alpha + beta + nj) + lgamma(alpha+yj) + lgamma(beta+nj-yj)
  interm <- sum(temp)+J*(lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta))
  pro <- exp(interm)
  pro*alpha*beta*(alpha+beta)^(-5/2)
}

lkrc <- function(nj, yj, alpha=1, beta=1) {
  J <- length(nj)
  if (J != length(yj)) stop("bad sizes")
  temp <- - lgamma(alpha + beta + nj) + lgamma(alpha+yj) + lgamma(beta+nj-yj)
  interm <- sum(temp)+J*(lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta))
  interm + log(alpha) + log(beta) + log(alpha+beta)*(-5/2)
}

# br is the object describing the posterior distribution
#
# get the expected values of the parameters in transformed x-y-space
expext <- function(br) {
  res <- br@results/sum(br@results)
  sx <- apply(res, 1, sum)
  ex <- sum(br@x*sx)
  sy <- apply(res, 2, sum)
  ey <- sum(br@y*sy)
  list(x=ex, y=ey)
}
# transform the parameters into interpretable alpha-beta-space.
xform <- function(xy) {
  x <- xy$x
  y <- xy$y
  alpha <- exp(x+y)/(1+exp(x))
  beta <- exp(y)/(1+exp(x))
  list(alpha=alpha, beta=beta, mean=alpha/(alpha+beta), size=alpha+beta)
}

##############################
# BetaRates is the main class, as well as the name of a constructor.

setClass("BetaRates",
         slots = c(k = 'numeric',
                   n = 'numeric',
                   x = 'numeric',
                   y = 'numeric',
                   results = 'matrix',
                   logresults = 'matrix'))

BetaRates <- function(k, n, x=seq(-3,3, length=100), y=x) {
  NX <- length(x)
  NY <- length(y)
  logresults <- matrix(NA, NX, NY)
  for (i in 1:NX) {
    for (j in 1:NY) {
      a <- exp( x[i] + y[j]) / (1 + exp(x[i]) )
      b <- exp( y[j]) / (1 + exp(x[i]) )
      logresults[i,j] <- lkrc(n, k, a, b)
    }
  }
  results <- exp(logresults - max(logresults))
  weight <- sum(results)
  results <- results/weight
  new("BetaRates", k=k, n=n, x=x, y=y, results=results, logresults=logresults)
}

##############################
# S4 methods
setMethod('summary', 'BetaRates', function(object, ...) {
  x <- xform(ee <- expext(object))
  cat("A BetaRates object constructed on data from", length(object@n),
      "different groups\n",
      "Posterior estimates of the parameters are:\n")
  print(unlist(ee))
  cat("These are equivalent to:\n")
  unlist(x)
})

# show the 2D posterior distribution on the transformed parameter space
setMethod("image", "BetaRates", function(x, col=greyscale(128), ...) {
  ii <- which(apply(x@results, 1, max)==max(x@results))
  jj <- which(apply(x@results, 2, max)==max(x@results))
  image(x@x, x@y, x@results, col=col, #$
        xlab="log(alpha/beta)", ylab="log(alpha+beta)", ...)
  points(x@x[ii], x@y[jj], pch=16, col='red')
  invisible(list(x=x@x[ii], y=x@y[jj]))
})


##############################
# TODO: Convert the retun value into its own class with
# some user-friendly methods
samplePosteriorRates <- function(br, nsamp=2000) {
  x <- br@x
  mu <- exp(x)/(1+exp(x))
  distri <- apply(br@results, 1, sum)/sum(br@results)

  temp <- runif(nsamp)
  xsamp <- unlist(sapply(temp, function(tt) {
    x[1+sum(cumsum(distri)<tt)]
  }))
  ysamp <- unlist(lapply(xsamp, function(i) {
    xi <- which(x==i)
    condist <- br@results[xi,]/sum(br@results[xi,])
    br@y[1+sum(cumsum(condist)<runif(1))] 
  }))
  beta <- exp(ysamp)/(1+exp(xsamp))
  alpha <- beta*exp(xsamp)
  theta <- matrix(NA, nrow=nsamp, ncol=length(br@n))
  for (i in 1:length(br@n)) {
    theta[,i] <- rbeta(nsamp, alpha+br@k[i], beta+br@n[i]-br@k[i])
  }
  theta <- data.frame(theta)
  brand <- names(br@n)
  if (!is.null(brand)) {
    if(all(!is.na(brand))) {
      colnames(theta) <- names(br@n)
    }
  }
  list(xy=data.frame(x=xsamp, y=ysamp), theta=theta)
}

##############################
# This is useful for deciding where to center the
# grid for computing the posterior distribution.
#
# Or better yet, figure out how to automatically set the grid.
guessCenter <- function(v) {
    m <- mean(v)
    s2 <- var(v)
    temp <- m*(1-m)/s2 - 1
    a <- m*temp
    b <- (1-m)*temp
    X <- log(a/b)
    Y <- log(a+b)
    list(X=X, Y=Y, a=a, b=b)
}

setClass("Deltas",
         slots = c(left = "numeric",
                   right = "numeric",
                   bottom = "numeric",
                   top = "numeric"))
setValidity("Deltas", function(object) {
  all( object@psi >= 0 ) & sum(object@psi) == 1
})
setMethod("initialize", "Deltas", function(.Object, left = 1, right = NULL,
                                           top = NULL, bottom = NULL, ...) {
  .Object <- callNextMethod(.Object) # in case this gets inherited
  inputs <- c(left, right, top, bottom)
  if (length(inputs) == 0) {
    inputs <- c(1,1,1,1)
  }
  while (length(inputs) < 4) {
    inputs <- c(inputs, inputs) # reuse earlier entries
  }
  if (length(inputs) != 4) { # only happens with 3 input values?
    stop("Incorrect number of delta values.")
  }
  .Object$left   <- inputs[1]
  .Object$right  <- inputs[2]
  .Object$top    <- inputs[3]
  .Object$bottom <- inputs[4]
  .Object
})
Deltas <- function(phi) {
  new("Deltas", psi = phi)
}
setAs("Deltas", "list", function(from) {
  list(left = from@left, right = from@right,
       top = from@top, bottom = from@bottom)
})

checkDeltas <- function(...) {
  as(Deltas(...), "list")
}

refineGrid <- function(K, N, x=seq(-3,3, length=100), y=x, deltas = NULL) {
  if (is.null(deltas)) {
    deltas <- list(left=1, right=1, bottom=1, top=1)
  }
  # can we really use an S4 class _just_ for error checking?
  deltas <- checkDeltas(deltas)
  guess <- guessCenter(K/N)
  repeat {
    br <- BetaRates(K, N, 
                    x = seq(guess$X - deltas$left, guess$X + deltas$right, length=100), 
                    y = seq(guess$Y - deltas$bottom, guess$Y + deltas$top, length=100))
    brr <- br@results
    edges <- 100*c(left=sum(brr[1,]), right=sum(brr[100,]), 
                   bottom=sum(brr[,1]), top=sum(brr[,100]))
    if (all(edges < 0.5)) break # no edge should have more than 0.5% of the mass
    temp <- names(which(edges >= 0.5))
    for (n in temp) deltas[[n]] <- deltas[[n]] + 0.5 # make a wider box
  }
  br
}
