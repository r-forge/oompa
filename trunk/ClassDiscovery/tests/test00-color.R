library(ClassDiscovery)
data <- matrix(1:1024, nrow=1024)
opar <- par(mfrow=c(2,3))
image(data, col=bluescale(64))
image(data, col=redscale(128))
image(data, col=greenscale(256))
image(data, col=redgreen(32))
image(data, col=blueyellow(64))
image(data, col=jetColors(64))
par(opar)
rm(data, opar) # cleanup

