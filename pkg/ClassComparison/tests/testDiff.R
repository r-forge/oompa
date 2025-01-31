library(ClassComparison)
sp <- options(scipen = 10)
nGenes <- 1000
nSamplesPerGroup <- 10
nGroups <- 2

set.seed(944637)
data <- matrix(rnorm(nGenes*nSamplesPerGroup*nGroups),
               nrow=nGenes)
classes <- factor(rep(c("A", "B"), each=nSamplesPerGroup))

mtt <- MultiTtest(data, classes)
summary(mtt, digits = 4)

suppressWarnings(mw <- MultiWilcoxonTest(data, classes))
summary(mw)

mlm <- MultiLinearModel(Y ~ classes, data.frame(classes=classes), data)
summary(mlm, digits = 4)

dud <- Dudoit(data, classes, nPerm=100, verbose=FALSE)
summary(dud, digits = 4)

tn <- TNoM(data, classes)
summary(tn, digits = 4)

sam <- Sam(data, classes)
summary(sam)

tgs <- TwoGroupStats(data, classes)
summary(tgs, digits = 4)
smoo <- SmoothTtest(tgs)
summary(smoo, digits = 4)

options(sp)
