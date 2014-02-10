# Copyright (C) Kevin R. Coombes, 2007-2012

## mosaic.R

setOldClass("hclust")

setClass('Mosaic',
         representation(data = 'matrix',
                        samples = 'hclust',
                        genes = 'hclust',
                        sampleMetric = 'character',
                        sampleLinkage = 'character',
                        geneMetric = 'character',
                        geneLinkage = 'character',
                        call = 'call',
                        name = 'character'))

Mosaic <- function(data, sampleMetric='pearson', sampleLinkage = 'average',
                   geneMetric='euclid', geneLinkage = 'average',
                   usecor=FALSE, center=FALSE, name='My mosaic') {
  call <- match.call()
  if(inherits(data, 'ExpressionSet')) {
    data <- exprs(data)
  }
  data <- as.matrix(as.data.frame(data))
  if(center == TRUE) {
    geneMean <- apply(data, 1, mean)
    data <- sweep(data, 1, geneMean, "-")
  }
  if(usecor == TRUE) {
    if(max(dim(data)) > 2000) {
      geneSD <- sqrt(apply(data, 1, var))
    }
    else {
      num.sample <- dim(data)[2]
      geneSD <- sqrt(diag(data %*% t(data)/(num.sample - 1)))
    }
    geneSD[geneSD <= 0] <- 1
    data <- sweep(data, 1, geneSD, "/")
  }
  samples <- hclust(distanceMatrix(data, sampleMetric), sampleLinkage)
  genes <- hclust(distanceMatrix(t(data), geneMetric), geneLinkage)
  val <- new('Mosaic', data=data, samples=samples, genes=genes,
             sampleMetric=sampleMetric, geneMetric=geneMetric,
             sampleLinkage=sampleLinkage, geneLinkage=geneLinkage,
             call=call, name=name)
}

setMethod('summary', signature(object='Mosaic'),
          function(object, ...) {
  cat(paste(object@name, ', an object of the ', class(object), ' class.\n\n', sep=''))
  cat(paste('Call:\n\t', as.character(list(object@call)), sep=''))
  cat(paste('\n\nSample dendrogram constructed with "', object@sampleLinkage,
            '" linkage and "', object@sampleMetric, '" distance metric.', sep=''))
  cat(paste('\n\nGene dendrogram constructed with "', object@geneLinkage,
            '" linkage and "', object@geneMetric, '" distance metric.\n', sep=''))
})

setMethod('plot', signature('Mosaic', 'missing'),
          function(x,
                   main=x@name,
                   center=FALSE,
                   scale=c('none', 'row', 'column'),
                   limits=NULL,
                   sampleColors=NULL,
                   sampleClasses=NULL,
                   geneColors=NULL,
                   geneClasses=NULL,
                   ...) {
  data <- x@data
  if(center == TRUE) {
    geneMean <- apply(data, 1, mean)
    data <- sweep(data, 1, geneMean, "-")
  }
  scale <- match.arg(scale)
  if (scale == 'row') {
    data <- t(scale(t(data)))
  } else if (scale == 'column') {
    data <- scale(data)
  }
  if (!is.null(limits)) {
    if (length(limits) < 2) {
      limits <- c(-limits, limits)
    }
    data[data > max(limits)] <- max(limits)
    data[data < min(limits)] <- min(limits)
  }
  arglist <- list(data, Colv=as.dendrogram(x@samples),
                  Rowv=as.dendrogram(x@genes), ...)
  if (!is.null(sampleClasses)) {
    if (length(sampleClasses) > 1) {
      classes <- as.numeric(sampleClasses)
    } else {
      classes <- cutree(x@samples, k=sampleClasses)
    }
    if (is.null(sampleColors)) {list
      sampleColors <- sample(colors(), max(classes))
    }
    arglist <- append(arglist, list(ColSideColors=sampleColors[classes]))
  }
  if (!is.null(geneClasses)) {
    if (length(geneClasses) > 1) {
      classes <- as.numeric(geneClasses)
    } else {
      classes <- cutree(x@genes, k=geneClasses)
    }
    if(min(classes)==0) classes <- 1 + classes
    if (is.null(geneColors)) {
      geneColors <- sample(colors(), max(classes))
    }
    arglist <- append(arglist, list(RowSideColors=geneColors[classes]))
  }
  arglist <- append(arglist, list(scale='none'))
  do.call('aspectHeatmap', arglist)
  invisible(arglist)
})

setMethod('pltree', signature(x='Mosaic'),
          function(x, colors, labels, ...) {
  if (missing(labels)) {
    labels <- dimnames(x@data)[[2]]
    if (is.null(labels)) {
      labels <- paste('X', 1:ncol(x@data), sep='')
    }
  }
  if (missing(colors)) {
    colors = rep('black', ncol(x@data))
  }
  plotColoredClusters(x@samples, labels, colors)
  invisible(x)
})

