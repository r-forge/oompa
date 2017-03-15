library(BetaModels)

# simulate some data
event <- c( 37,  4,  6,  1,  2, 10,  1, 13,   7,  1,  10)
total <- c(137, 18, 18, 26, 24, 45, 12, 43, 162, 78, 280)
bg <- total - event

# Use a chi-squared test for overall significance
chisq.test(data.frame(event, bg))
chisq.test(data.frame(event, bg), simulate.p.value=TRUE)

# Haven't yet figured out to automate this step. But not weilling
# to commit to making this method part of the documented interface.
BetaModels:::guessCenter(event/total)

# compute the posterior
br <- BetaRates(event, total, x=seq(-3, 0, length=100), y=seq(0, 3, length=100))
# view it
image(br)
# summarize it
summary(br)

# sample from the posterior distribution
spr <- samplePosteriorRates(br, nsamp=2000)
theta <- spr$theta

# observed point estimate of rates
round(naive <- event/total, 4)
# bayesian mean rate
round(bayes <- apply(theta, 2, mean), 4)
# visualize "shrinking toward the mean"
plot(naive, bayes, pch=16, xlim=c(0, 0.35), ylim=c(0, 0.35))
abline(0,1)
abline(h=sum(event)/sum(total), col='blue')
abline(h=mean(event/total), col='purple')

# visualize the posterior distribution
boxplot(theta)

if(require(vioplot)) {
  vioplot(theta[,1], theta[,2], theta[,3], theta[,4],
          theta[,5], theta[,6], theta[,7], theta[,8],
          theta[,9], theta[,10], theta[,11], col='skyblue'
          )
}

# posterior probability that rate in group 3 is higher than other groups
j <- 3
sapply(1:11, function(i) mean(theta[,j] > theta[,i]))

# posterior probability that rate in group 3 is at least 5% higher than other groups
j <- 3
sapply(1:11, function(i) mean(theta[,j] > theta[,i] + 0.05))

# marginal posterior distributions
hist(spr$xy$x, breaks=55)
hist(spr$xy$y, breaks=55)

ab <- as.data.frame(BetaModels:::xform(spr$xy))
hist(ab$alpha, breaks=55)
hist(ab$beta, breaks=55)
hist(ab$mean, breaks=55)
hist(ab$size, breaks=55)
