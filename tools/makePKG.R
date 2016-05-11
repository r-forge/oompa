library(tools)

thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
choppedRVer <- gsub("(\\w+).(\\w+).(\\w+)", "\\1.\\2", thisRVer)

home <- getwd()
setwd(paste("Binary", thisRVer, "x64", sep='-'))
write_PACKAGES(".", type="win.binary")

setwd(home)
setwd(paste("Build", thisRVer, sep='-'))
write_PACKAGES(".", type="source")
