library(BetaModels)

set.seed(73892)
datavec <- c(rbeta(100, 1, 1),
             rbeta(200, 7, 4))

# randomly initialize Z
temp <- sample(2, length(datavec), replace = TRUE)
Z <- matrix(0, nrow = length(datavec), ncol = 2)
for (I in 1:nrow(Z)) Z[I, temp[I]] <- 1

# initialize parameters (2 identical components)
initparam <- rep(1/2, 4)
NegBetaLogLike(initparam, datavec, Z)

# use true parameters
NegBetaLogLike(c(1, 1, 7, 4), datavec, Z)


model <- BetaMixture(datavec, debug = TRUE, forever = 100)
summary(model)
hist(model, breaks=35)

