library(BetaModels)

set.seed(73892)
datavec <- c(rbeta(130, 1, 4),
             rbeta(170, 7, 4))

# randomly initialize Z
temp <- sample(2, 200, replace = TRUE)
Z <- matrix(0, nrow = 100, ncol = 2)
for (I in 1:100) Z[I, temp[I]] <- 1

# initialize parameters (2 identical components)
initparam <- rep(1/2, 4)
NegBetaLogLike(initparam, datavec, Z)

# use true parameters
NegBetaLogLike(c(1, 4, 7, 4), datavec, Z)


model <- BetaMixture(datavec, debug = TRUE, forever = 10)
summary(model)
hist(model, breaks=35)

