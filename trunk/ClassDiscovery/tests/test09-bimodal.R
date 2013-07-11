library(ClassDiscovery)
# want to test the bimodal index method
set.seed(632984) # for reproducibility
# start by generating a dataset with no bimodality
nGenes <- 500
nSamples <- 200
temp <- matrix(rnorm(nGenes*nSamples), nrow=nGenes)
# define the number of samples in the smaller group
fracs <- rep(nSamples*seq(0.1, 0.5, by=0.1), each=20)
# define the spread between the two means, per gene
delta <- rnorm(length(fracs), 3,1)
# now shift the means for the smaller group. This is
# woefully inefficient.
for (i in 1:length(fracs)) {
  temp[i, 1:fracs[i]] <- temp[i, 1:fracs[i]] + delta[i]
}
# now we compute the bimodal index
bim <- bimodalIndex(temp)
summary(bim)
summary(bim[1:100,])
head(bim)
round(bim[1:100,], digits=6)
