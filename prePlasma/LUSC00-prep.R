## Load the ESCA TCGA data
clean  <- "F:/TCGA/clean"
intern <- list(Clinical = "clin",
               CNV = "segments",
               MAF = c("mutCall", "mutations", "mutFraction"),
               Meth450 = "stuff",
               miRSeq = c("lcounts", "ncounts"),
               mRNASeq = c("lcounts", "ncounts"),
               RPPA = "ncounts")
cancer  <- "LUSC"

for(TY in names(intern)) {
  cat(TY,"\n", file = stderr())
  load(file.path(clean, TY, paste(cancer, "Rda", sep = ".")))
  assign(TY, get(intern[[TY]][1]), .GlobalEnv)
  rm(list = intern[[TY]])
}
rm(cancer, clean, TY, intern)

## Separate clinical covariates from outcomes.
temp <- Clinical
foo <- toupper(temp$patient_id)
rownames(temp) <- foo
Outcome <- temp[, 2:5]
Clinical <- temp[, -(1:5)]

## Mak sense of the Outcome data
summary(Outcome)
time <- as.integer(as.character(Outcome$days_to_death))
Outcome$days_to_death <- time
lfu <- as.integer(as.character(Outcome$days_to_last_followup))
Outcome$days_to_last_followup <- lfu
time[is.na(time)] <- lfu[is.na(time)]
Outcome$Days <- time
summary(Outcome)
## Need to get rid of samples where we don't know the outcome
Clinical <- Clinical[!is.na(Outcome$Days),]
Outcome <- Outcome[!is.na(Outcome$Days),]
summary(Outcome)

## Remove things with lots of NA's
nacount <- apply(Clinical, 2, function(x) sum(is.na(x)))
plot(sort(nacount))
abline(h = 126) # 25%
colnames(temp)[nacount > 126]
colnames(temp)[nacount < 126]
Clinical <- Clinical[, nacount < 126]

## Remove things with only one unique non-trivial value
nvals <- apply(Clinical, 2, function(X) {
  X <- X[!is.na(X)]
  length(unique(X))
})
Clinical <- Clinical[, nvals > 1]
nvals <- apply(Clinical, 2, function(X) {
  X <- X[!is.na(X)]
  length(unique(X))
})
table(nvals)
rm(nvals)

## Remove redundant columns, and columns with useless lab information
summary(Clinical)
numerize <- c(1, 16)
factors <- c(2, 5, 9:11, 12, 17, 20, 22:25, 27:28, 30:31, 37:39)
dups <- c(7, 13, 18:19, 21, 29, 34:36, 40:46)
throwaway <- c(3:4, 6, 8, 14:15, 26, 32:33)

compare <- function(A, B) {
  table(Clinical[, A], Clinical[, B])
}
compare(10, 13)
compare(28, 29)
compare(31, 2)
compare(34, 1)
compare(20, 36)
compare(21, 37)
compare(19, 38)
compare(18, 39)
compare(9, 40)

colnames(Clinical)[numerize]
colnames(Clinical)[factors]
colnames(Clinical)[dups]
colnames(Clinical)[throwaway]

Clinical <- Clinical[, c(numerize, factors)]
for (I in 1:2) {
  Clinical[,I] <- as.numeric(as.character(Clinical[,I]))
}

summary(Clinical)
rm(temp, foo, nacount, I, compare, numerize, factors, dups, throwaway)

###########################################
## Clean up MAF data
dim(MAF)
X <- t(MAF)
sum(is.na(X))
sum(X == "")
Y <- 1*(X != "WT")
Y[X== ""] <- NA
foo <- sapply(strsplit(rownames(Y), "-"), function(W) W[3])
goo <- substring(sapply(strsplit(rownames(Y), "-"), function(W) W[4]), 1, 2)
Y <- Y[goo == "01",]
rownames(Y) <- foo
MAF <- t(Y)
knock <- apply(MAF, 1, function(X) sum(is.na(X)))
table(knock)
mu <- apply(MAF, 1, mean, na.rm=TRUE)
M0 <- 0.04
sum(mu > M0)
MAF <- MAF[mu > M0,]
rm(X, Y, foo, mu, knock)
MAF <- MAF[, colnames(MAF) %in% rownames(Outcome)]
dim(MAF)

## Clean up data from Methylation 450 arrays
dim(Meth450$betas)
m450info <- Meth450$info
m450beta <- Meth450$betas
goo <- substring(sapply(strsplit(colnames(m450beta), "\\."), function(W) W[4]), 1, 2)
table(goo)
m450beta <- m450beta[, goo %in% c("01", "06")]
#m450beta <- m450beta[, -156]
foo <- sapply(strsplit(colnames(m450beta), "\\."), function(W) W[3])
colnames(m450beta) <- foo
rm(foo)
mu <- apply(m450beta, 1, mean, na.rm = TRUE)
sigma <- apply(m450beta, 1, sd, na.rm = TRUE)
M0 <- 0.15
S0 <- 0.25
smoothScatter(mu, sigma)
abline(v = c(M0, 1 - M0), h = S0, col = "orange", lwd=2)
keep <- !is.na(mu) & mu > M0 & mu < (1-M0) & sigma > S0
summary(keep)
m450info <- m450info[keep,]
m450beta <- m450beta[keep,]
m450beta <- m450beta[, colnames(m450beta) %in% rownames(Outcome)]
dim(m450beta)

## Clean up miRSeq data
dim(miRSeq)
goo <- substring(sapply(strsplit(colnames(miRSeq), "\\."), function(W) W[4]), 1, 2)
miRSeq <- miRSeq[, goo %in% c("01", "06")]
foo <- sapply(strsplit(colnames(miRSeq), "\\."), function(W) W[3])
colnames(miRSeq) <- foo
rm(foo)
range(miRSeq)
mu <- apply(miRSeq, 1, mean, na.rm = TRUE)
sigma <- apply(miRSeq, 1, sd, na.rm = TRUE)
S0 <- 0.1
smoothScatter(mu, sigma)
abline(h = S0)
keep <- sigma > S0
summary(keep)
miRSeq <- miRSeq[keep,]
rm(mu, sigma, keep)
miRSeq <- miRSeq[, colnames(miRSeq) %in% rownames(Outcome)]
dim(miRSeq)

## Clean up mRNA data
dim(mRNASeq)
goo <- substring(sapply(strsplit(colnames(mRNASeq), "\\."), function(W) W[4]), 1, 2)
table(goo)
mRNASeq <- mRNASeq[, goo %in% c("01", "06")]
foo <- sapply(strsplit(colnames(mRNASeq), "\\."), function(W) W[3])
colnames(mRNASeq) <- foo
rm(foo)
range(mRNASeq)
mu <- apply(mRNASeq, 1, mean, na.rm = TRUE)
sigma <- apply(mRNASeq, 1, sd, na.rm = TRUE)
M0 <- 5
S0 <- 1.25
smoothScatter(mu, sigma)
abline(v=M0, h=S0, col="yellow", lwd=2)
keep <- mu > M0 & sigma > S0
summary(keep)
mRNASeq <- mRNASeq[keep,]
rm(mu, sigma, keep)
mRNASeq <- mRNASeq[, colnames(mRNASeq) %in% rownames(Outcome)]
dim(mRNASeq)

## Clean up RPPA data
dim(RPPA)
goo <- substring(sapply(strsplit(colnames(RPPA), "\\."), function(W) W[4]), 1, 2)
table(goo)
RPPA <- RPPA[, goo %in% c("01", "06")]
foo <- sapply(strsplit(colnames(RPPA), "\\."), function(W) W[3])
colnames(RPPA) <- foo
rm(foo)
range(RPPA)
mu <- apply(RPPA, 1, mean, na.rm = TRUE)
sigma <- apply(RPPA, 1, sd, na.rm = TRUE)
smoothScatter(mu, sigma)
RPPA <- RPPA[, colnames(RPPA) %in% rownames(Outcome)]
rm(goo, mu, sigma)
dim(RPPA)

## Assemble all the data sets
assemble <- list(Clinical = t(Clinical),
                 MAF = MAF,
                 Meth450 = m450beta,
                 miRSeq = miRSeq,
                 mRNASeq = mRNASeq,
                 RPPA = RPPA)
sapply(assemble, dim)


save(Outcome, assemble, m450info, file = "TCGA-LUSC0.RData")

## verify that it can be re-loaded
rm(assemble)
load("TCGA-LUSC0.RData")
