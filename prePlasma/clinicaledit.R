library(plasma)
load("TCGA-ESCA0.Rdata")
clin <- assemble$Clinical
clindf  <-  as.data.frame(t(clin))

## separate into continuous and categorical
g <- apply(clindf, 2, as.numeric)
idx_cat <- which(apply(g, 2, function(x) sum(is.na(x)) == nrow(g)))
idx_num <- which(apply(g, 2, function(x) sum(is.na(x)) != nrow(g)))
clincat <- clindf[,idx_cat]
clinnum <- clindf[,idx_num]

## create continuous clinical variables
clincont <- apply(clinnum, 2, as.numeric)
rownames(clincont) <- rownames(clinnum)
## peek at it
dim(clincont)
summary(clincont)
## one strongly suspects that "tobacco_smoking_history" is categorical, not continuous
## But it is critical to know what the numbers mean. This is from
##    https://github.com/PoisonAlien/maftools/issues/410
toblevels <- c("NeverSmoker", "CurrentSmoker", "FormerSmokerGreater15",
               "FormerSmokerLess15", "FormerSmokerUnknown", NA, NA)
w <- which(colnames(clincont) == "tobacco_smoking_history")
tobaccofactor <- factor(toblevels[clincont[,w]], levels = toblevels[1:4])
clincont <- clincont[, -w]
rm(w)
summary(clincont)
clincat$tobacco <- tobaccofactor

## create binary clinical variables
## KRC: How do we keep track of what "1" means after these transformations?
##      It appears to be encoded as part of the new expanded column names
makeBinary <- function(x){
  options(na.action = 'na.pass')
  BIN <- factor(x)
  b <- model.matrix(~ 0 + BIN)
  colnames(b)  <-  gsub("BIN", "", colnames(b))
  return(b)
  options(na.action = 'na.omit')
}
clincat1 <- as.data.frame(apply(clincat, 2, function(x) makeBinary(x)))
## get rid of all the columns with "no" in it since it's redundant from the "yes" columns
clincat2 <- clincat1[ !(colnames(clincat1) %in% grep(".no", colnames(clincat1), value=TRUE))]
## drop the ".yes" in the remaining columns
colnames(clincat2) <- gsub(".yes", "", colnames(clincat2))
rownames(clincat2) <- rownames(clincat)
## peek
pos <- apply(clincat2, 2, function(X) sum(X == 1, na.rm = TRUE))
neg <- apply(clincat2, 2, function(X) sum(X == 0, na.rm = TRUE))
cbind(pos, neg)
## 1. Do we really need to track uk and usa vversiosn of barrets esophagus? Apparently, yes;
##    since the diagnostic criteria are different.
## 2. Do we want to track city where the samples was accrued? What about country?
## 3. What about all those variants of stages or partial stages (M, N, T)?
## Let's get rid of the cities
g <- grepl("city", colnames(clincat2))
clincat3 <- clincat2[, !g]
dim(clincat3)
## and let's keep "pathology" but mnot "pathologic" stage components
g <- grepl("pathologic", colnames(clincat3))
clincat4 <- clincat3[, !g]
dim(clincat4)
## peek
pos <- apply(clincat4, 2, function(X) sum(X == 1, na.rm = TRUE))
neg <- apply(clincat4, 2, function(X) sum(X == 0, na.rm = TRUE))
cbind(pos, neg)
## Let's also get rid of country, since so many levels only appear once
g <- grepl("country", colnames(clincat4))
clincat5 <- clincat4[, !g]
dim(clincat5)
## peek
pos <- apply(clincat5, 2, function(X) sum(X == 1, na.rm = TRUE))
neg <- apply(clincat5, 2, function(X) sum(X == 0, na.rm = TRUE))
cbind(pos, neg)



clinbin  <-  t(clincat5)
clincont <- t(clincont)

assemble <- list(ClinicalBin = clinbin,
                 ClinicalCont = clincont,
                 MAF = assemble$MAF,
                 Meth450 = assemble$Meth450,
                 miRSeq = assemble$miRSeq,
                 mRNASeq = assemble$mRNASeq,
                 RPPA = assemble$RPPA)
sapply(assemble, dim)
save(Outcome, assemble, m450info, file = "TCGA-ESCA-full.RData")

set.seed(97531)
miRSeq <- assemble$miRSeq[, sample(ncol(assemble$miRSeq), 166)]
assemble$miRSeq <- miRSeq
set.seed(24680)
mRNASeq <- assemble$mRNASeq[, sample(ncol(assemble$mRNASeq), 157)]
assemble$mRNASeq <- mRNASeq
sapply(assemble, dim)

save(Outcome, assemble, m450info, file = "TCGA-ESCA1.RData")
