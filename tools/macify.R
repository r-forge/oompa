library(tools)
library(utils)

thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
choppedRVer <- gsub("(\\w+).(\\w+).(\\w+)", "\\1.\\2", thisRVer)

home <- getwd()
setwd("../..")
if (!file.exists("Working")) dir.create("Working")
setwd("Working")

winbin <- file.path("..", paste("Binary", thisRVer, "x64", sep='-'))
macbin <- file.path("..", paste("MacBinary", thisRVer, "x64", sep='-'))

if (!file.exists(winbin)) stop("no windows binarues")
if (!file.exists(macbin)) dir.create(macbin)

files <- dir(winbin, pattern="zip")
for (myfile in files) {
  cat("xforming", myfile, "\n", file=stderr())
  unzip(file.path(winbin, myfile))
  base <- strsplit(myfile, "_")[[1]][1]
  if (!file.exists(base)) stop("Cannot find directory", base)
  tarfile <- file.path(macbin, sub("zip", "tgz", myfile))
  tar(tarfile, base, compression="gzip")
  unlink(base, recursive=TRUE)
}
