# Copyright (C) Kevin R. Coombes, 2007-2017

redscale <- function(N) {
  rgb((1:N)-1, 0, 0, maxColorValue=N)
}

greenscale <- function(N) {
  rgb(0, (1:N)-1, 0, maxColorValue=N)
}

bluescale <- function(N) {
  rgb(0, 0, (1:N)-1, maxColorValue=N)
}

blueyellow <- function(N) {
  x <- (1:N)-1
  rgb(x, x, rev(x), maxColorValue=N)
}

cyanyellow <- function(N) {
  x <- colorRampPalette(c("cyan", "black", "yellow"))
  x(N)
}

greyscale <- function (N) 
{
  x <- (1:N)-1
  rgb(x, x, x, maxColorValue = N)
}
grayscale <- greyscale

## Fix to redgreen contributed by Karl Kashofer
## Fixed yet again, since it got the maxcolorvalue wrong
## and thus washed everything out.
redgreen <- function (N) {
  if (N%%2 == 0) {
    A  <- N/2
    r <- c(rep(0, A), 0:A)[-A]
    g <- -c(-A:0, rep(0, A))[-A-1]
  } else {
    A  <- (N-1)/2
    r <- c(rep(0, A), 0:A)
    g <- -c(-A:0, rep(0, A))
  }
  rgb(r, g, rep(0, N), maxColorValue = A)
}


## jet Colors from Keith Baggerly, based on MATLAB default color map
jetColors <- function(N){
  k <- ceiling(N/4)
  temp.red <- c(rep(0,2*k), 1:k, rep(k,k-1), k:1)
  temp.green <- c(rep(0,k), 1:k, rep(k,k-1), k:1, rep(0,k))
  temp.blue <- c(1:k, rep(k,k-1), k:1, rep(0,2*k))
  temp.rgb <- cbind(temp.red, temp.green, temp.blue)
  delta <- 5*k-1 - N
  delta <- ceiling(delta/2)
  temp.rgb <- temp.rgb[delta:(delta+N-1),]/k

  ## assemble everything last value is returned
  rgb(temp.rgb[,1], temp.rgb[,2], temp.rgb[,3])

}

################################################################
# some common colors that we use for specialized plots

oompaColor = list(
  BORING              = 'gray',
  SIGNIFICANT         = 'red',
  EXPECTED            = 'blue',
  OBSERVED            = 'darkgreen',
  
  PERMTEST            = 'magenta',
  FITTED              = 'orange',
  CENTRAL.LINE        = 'blue',
  CONFIDENCE.CURVE    = 'red3',

  BAD.REPLICATE       = 'purple1',
  WORST.REPLICATE     = 'purple3',
  FOLD.DIFFERENCE     = 'skyblue',
  BAD.REPLICATE.RATIO = 'violetred',

  TOP.TEN             = 'cadetblue',
  BOTTOM.TEN          = 'pink',
  TOP.TEN.SOLO        = 'palegreen',
  BOTTOM.TEN.SOLO     = 'deeppink'
  )

