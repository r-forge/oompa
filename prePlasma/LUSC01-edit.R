library(plasma)
load("TCGA-LUSC0.Rdata")
clin <- assemble$Clinical
clindf  <-  as.data.frame(t(clin))
w <- which(colnames(clindf) == "tissue_source_site")
clindf <- clindf[, -w]

## separate into continuous and categorical
g <- apply(clindf, 2, as.numeric)
idx_cat <- which(apply(g, 2, function(x) sum(is.na(x)) == nrow(g)))
idx_num <- which(apply(g, 2, function(x) sum(is.na(x)) != nrow(g)))
idx_num <- 1:2
idx_cat <- 3:ncol(clindf)
clincat <- clindf[,idx_cat]
clinnum <- clindf[,idx_num]

## create continuous clinical variables
clincont <- apply(clinnum, 2, as.numeric)
rownames(clincont) <- rownames(clinnum)
## peek at it
dim(clincont)
summary(clincont)

## make everyt5hing factors
for (I in 1:ncol(clincat)) {
  clincat[, I] <- as.factor(clincat[, I])
}
temp <- as.character(clincat$history_of_neoadjuvant_treatment)
temp[temp != "no"] <- "yes"
clincat$history_of_neoadjuvant_treatment<- as.factor(temp)
temp <- as.character(clincat$other_dx)
temp[!is.na(temp) & temp != "no"] <- "yes"
clincat$other_dx<- as.factor(temp)
temp <- as.character(clincat$race)
temp[!is.na(temp) & temp != "white"] <- "other"
clincat$race<- as.factor(temp)
## one strongly suspects that "tobacco_smoking_history" is categorical, not continuous
## But it is critical to know what the numbers mean. This is from
##    https://github.com/PoisonAlien/maftools/issues/410
toblevels <- c("NeverSmoker", "CurrentSmoker", "FormerSmokerGreater15",
               "FormerSmokerLess15", "FormerSmokerUnknown", NA, NA)
w <- which(colnames(clincat) == "tobacco_smoking_history")
tobaccofactor <- factor(toblevels[clincat[,w]], levels = toblevels[1:5])
clincat[, w] <- tobaccofactor
rm(w, tobaccofactor)
summary(clincat)


## create binary clinical variables
## KRC: How do we keep track of what "1" means after these transformations?
##      It appears to be encoded as part of the new expanded column names
makeBinary <- function(x){
  options(na.action = 'na.pass')
  on.exit(options(na.action = 'na.omit'))
  BIN <- factor(x)
  b <- model.matrix(~ 0 + BIN)
  colnames(b)  <-  gsub("BIN", "", colnames(b))
  return(b)
}
clincat1 <- as.data.frame(apply(clincat, 2, function(x) makeBinary(x)))
## get rid of all the columns with "no" in it since it's redundant from the "yes" columns
clincat2 <- clincat1[ !(colnames(clincat1) %in% grep(".no", colnames(clincat1), value=TRUE))]
## drop the ".yes" in the remaining columns
colnames(clincat2) <- gsub(".yes", "", colnames(clincat2))
rownames(clincat2) <- rownames(clincat)
summary(clincat2)


clinbin  <-  t(clincat2)
clincont <- t(clincont)

assemble <- list(ClinicalBin = clinbin,
                 ClinicalCont = clincont,
                 MAF = assemble$MAF,
                 Meth450 = assemble$Meth450,
                 miRSeq = assemble$miRSeq,
                 mRNASeq = assemble$mRNASeq,
                 RPPA = assemble$RPPA)
sapply(assemble, dim)
save(Outcome, assemble, m450info, file = "TCGA-LUSC1.RData")
