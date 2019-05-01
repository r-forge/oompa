suppressWarnings( RNGversion("3.5.3") )
set.seed(932451)
#################################################################
# load the data set
library(ClassDiscovery)
require(oompaData) # for the prostate cancer data set with 2000 genes
data(expression.data)
data(clinical.info)
# make sure the Status variable has "N"ormals as the first level
# Note that "T" = primary prostate tumor and "L" = lymph node metastasis
temp <- ordered(clinical.info$Status, c('N', 'T', 'L'))
clinical.info$Status <- temp
rm(temp)

#################################################################
# create an exprSet
suppressMessages( require(Biobase) )
# create a phenoData object and an exprSet
vl <- data.frame(labelDescription=dimnames(clinical.info)[[2]])
rownames(vl) <- as.character(vl[,1])
pheno <- new('AnnotatedDataFrame',
             data=clinical.info,
             varMetadata=vl)
es <- new('ExpressionSet',
          phenoData=pheno,
          exprs=as.matrix(expression.data))
# we don't need the original data now that it has been
# incorporated into the exprSet
rm(expression.data)
rm(vl, pheno)

#################################################################
# Now we can start exercising the ClassDiscovery package

##################################
#windows(width=14, height=7, pointsize=10)
par(mfrow=c(2,2))
cluster3(es)
par(mfrow=c(1,1))

##################################
spc <- SamplePCA(es, 'Status')
levels(spc@splitter)
colorScheme <- c('green', 'magenta', 'red')
plot(spc, col=colorScheme)
legend(-25, -25, levels(spc@splitter), pch=15, col=colorScheme)

##################################
metric <- 'pearson'
linkage <- 'complete'
hc <- hclust(distanceMatrix(exprs(es), metric), linkage)
plotColoredClusters(hc, labs=colnames(exprs(es)),
                    col=colorScheme[as.numeric(spc@splitter)])

##################################
if (FALSE) { # don't run; this is slow
bc <- BootstrapClusterTest(es, cutHclust, k = 12, nTimes=200, verbose=FALSE,
                           metric=metric, method=linkage)
summary(bc)
hist(bc)
image(bc, dendrogram=hc, col=blueyellow(64))
}

##################################
if (FALSE) { # don't run;  this is slow
pc <- PerturbationClusterTest(es, cutHclust, k = 10, nTimes=100, verbose=FALSE,
                           noise=1, metric=metric, method=linkage)
summary(pc)
hist(pc)
image(pc, dendrogram=hc, col=blueyellow(64))
}

##################################
# mosaic

mos <- Mosaic(es,center=TRUE, usecor=TRUE, geneMetric='pearson')

dimnames(pData(es))[[2]]

plot(mos, sampleClasses=pData(es)[,'Status'], sampleColors=colorScheme,
     col=blueyellow(64), limits=c(-1,1), geneClasses=8)

if (FALSE) { # don't really need more heatmaps
plot(mos, sampleClasses=pData(es)[,'ChipType'], sampleColors=colorScheme,
     col=blueyellow(64), limits=c(-1,1), geneClasses=8)

plot(mos, sampleClasses=pData(es)[,'Subgroups'], sampleColors=colorScheme,
     col=blueyellow(64), limits=c(-1,1), geneClasses=8)
}

##################################
# clean everything up

rm(es, spc, colorScheme, hc, mos)
