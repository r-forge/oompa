library(plasma)
data("TCGA-ESCA")
clin<-assemble$Clinical
clindf<-as.data.frame(t(clin))

# separate into continuous and categorical
g<-apply(clindf, 2, as.numeric)
idx_cat<-which(apply(g, 2, function(x) sum(is.na(x))==nrow(g)))
idx_num<-which(apply(g, 2, function(x) sum(is.na(x))!=nrow(g)))
clincat<-clindf[,idx_cat]
clinnum<-clindf[,idx_num]

# create continuous clinical variables
clincont<-t(apply(clinnum, 2, as.numeric))
colnames(clincont)<-rownames(clinnum)

# create binary clinical variables
makeBinary<-function(x){
  options(na.action='na.pass')
  BIN = factor(x)
  b=model.matrix(~0+BIN)
  colnames(b) <- gsub("BIN","",colnames(b))
  return(b)
  options(na.action='na.omit')
}
clincat1<-as.data.frame(apply(clincat, 2, function(x) makeBinary(x)))
# get rid of all the columns with "no" in it since it's redundant from the "yes" columns
clincat2<-clincat1[ !(colnames(clincat1) %in% grep(".no", colnames(clincat1), value=TRUE))]
# drop the ".yes" in the remaining columns
colnames(clincat2)<-gsub(".yes", "", colnames(clincat2))
rownames(clincat2)<-rownames(clincat)
clinbin <- t(clincat2)