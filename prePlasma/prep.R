clean  <- "F:/TCGA/clean"
intern <- list(Clinical = "clin",
               CNV = "segments",
               MAF = c("mutCall", "mutations", "mutFraction"),
               Meth450 = "stuff",
               miRSeq = c("lcounts", "ncounts"),
               mRNASeq = c("lcounts", "ncounts"),
               RPPA = "ncounts")
cancer  <- "ESCA"

for(TY in names(intern)) {
  cat(TY,"\n", file = stderr())
  load(file.path(clean, TY, paste(cancer, "Rda", sep = ".")))
  assign(TY, get(intern[[TY]][1]), .GlobalEnv)
  rm(list = intern[[TY]])
}


temp <- Clinical
foo <- toupper(temp$patient_id)
rownames(temp) <- foo
Outcome <- temp[, 2:5]

nacount <- apply(Clinical, 2, function(x) sum(is.na(x)))
plot(sort(nacount))
abline(h=53)

Clinical <- temp[, nacount < 53]
Clinical <- Clinical[, c(4:6, 10:12, 15, 19:20, 22:25, 28,
                         32:33, 36, 39:48, 50:53, 58, 60, 71:73, 79)]
rm(temp, foo)


X <- t(MAF)
sum(is.na(X))
sum(X == "")
Y <- 1*(X != "WT")
Y[is.na(X)] <- NA
foo <- sapply(strsplit(rownames(Y), "-"), function(W) W[3])
goo <- substring(sapply(strsplit(rownames(Y), "-"), function(W) W[4]), 1, 2)
Y <- Y[goo == "01",]
rownames(Y) <- foo
MAF <- t(Y)
mu <- apply(MAF, 1, mean, na.rm=TRUE)
sum(mu > 0.01)
MAF <- MAF[mu > 0.01,]
rm(X, Y, foo, mu)
MAF <- MAF[, colnames(MAF) %in% rownames(Outcome)]

m450info <- Meth450$info
m450beta <- Meth450$betas
goo <- substring(sapply(strsplit(colnames(m450beta), "\\."), function(W) W[4]), 1, 2)
table(goo)
m450beta <- m450beta[, goo %in% c("01", "06")]
m450beta <- m450beta[, -156]
foo <- sapply(strsplit(colnames(m450beta), "\\."), function(W) W[3])
colnames(m450beta) <- foo
rm(foo)
mu <- apply(m450beta, 1, mean, na.rm = TRUE)
sigma <- apply(m450beta, 1, sd, na.rm = TRUE)
smoothScatter(mu, sigma)
abline(h = 0.02, v=0.1, col = "orange", lwd=2)
keep <- !is.na(mu) & mu > 0.02 & sigma > 0.1
summary(keep)
m450info <- m450info[keep,]
m450beta <- m450beta[keep,]
m450beta <- m450beta[, colnames(m450beta) %in% rownames(Outcome)]

goo <- substring(sapply(strsplit(colnames(miRSeq), "\\."), function(W) W[4]), 1, 2)
miRSeq <- miRSeq[, goo %in% c("01", "06")]
foo <- sapply(strsplit(colnames(miRSeq), "\\."), function(W) W[3])
colnames(miRSeq) <- foo
rm(foo)
range(miRSeq)
mu <- apply(miRSeq, 1, mean, na.rm = TRUE)
sigma <- apply(miRSeq, 1, sd, na.rm = TRUE)
smoothScatter(mu, sigma)
abline(h=0.05)
keep <- sigma > 0.05
summary(keep)
miRSeq <- miRSeq[keep,]
rm(mu, sigma, keep)
miRSeq <- miRSeq[, colnames(miRSeq) %in% rownames(Outcome)]

goo <- substring(sapply(strsplit(colnames(mRNASeq), "\\."), function(W) W[4]), 1, 2)
table(goo)
mRNASeq <- mRNASeq[, goo %in% c("01", "06")]
foo <- sapply(strsplit(colnames(mRNASeq), "\\."), function(W) W[3])
colnames(mRNASeq) <- foo
rm(foo)
range(mRNASeq)
mu <- apply(mRNASeq, 1, mean, na.rm = TRUE)
sigma <- apply(mRNASeq, 1, sd, na.rm = TRUE)
smoothScatter(mu, sigma)
abline(v=4, h=0.7, col="orange", lwd=2)
keep <- mu > 4 & sigma > 0.7
summary(keep)
mRNASeq <- mRNASeq[keep,]
rm(mu, sigma, keep)
mRNASeq <- mRNASeq[, colnames(mRNASeq) %in% rownames(Outcome)]

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

assemble <- list(Clinical = t(Clinical),
                 MAF = MAF,
                 Meth450 = m450beta,
                 miRSeq = miRSeq,
                 mRNASeq = mRNASeq,
                 RPPA = RPPA)
sapply(assemble, dim)
