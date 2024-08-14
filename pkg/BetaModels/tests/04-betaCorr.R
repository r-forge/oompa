library(BetaModels)

set.seed(12345)
cc <- c(rbeta(4600, 24, 24), rbeta(400, 8, 8))
rr <- 2*cc-1
fit <- ebCorrelation(rr, 51)
hist(fit)
plot(fit, prior = 0.85)
countSignificant(fit, prior = 0.85, significance = 0.8)
cutoffSignificant(fit, prior = 0.85, significance = 0.8)
summary(fit)
summary(fit, prior = 0.85, significance = 0.8)
