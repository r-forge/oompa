library(colorspace)
temp <- as.matrix(read.table("glasbey.txt", sep="\t", header=TRUE))
glasbey <- hex(sRGB(temp/255))
names(glasbey) <- paste("G", 1:32, sep='')
save(glasbey, file="glasbey.rda")
