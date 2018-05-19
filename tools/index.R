library(tools)
ver <- paste(R.version$major, R.version$minor, sep='.')
setwd("../..")

bin <- paste("Binary", ver, "x64", sep='-')
write_PACKAGES(bin, type="win.binary", verbose=TRUE)

src <- paste("Build", ver, sep='-')
write_PACKAGES(src, type="source", verbose=TRUE)

mac <- paste("BinaryMac", ver, "x64", sep='-')
write_PACKAGES(mac, type="mac.binary", verbose=TRUE)
