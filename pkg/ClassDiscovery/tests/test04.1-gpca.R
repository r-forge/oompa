library(ClassDiscovery)
RNGversion("3.5.3")
set.seed(341392)
# simulate samples from thre different groups, with generic genes
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)

# perform PCA in gene space
gpc <- GenePCA(dd)

# plot the results
plot(gpc)

# cleanup
rm(d1, d2, d3, dd, gpc)
