###
### oompaLite.R
###

oompaROOT <- "http://silicovore.com/OOMPA"

sourceOOMPAinstallScript <- function() {
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    ## Verify we're running a recent enough version of R.
    requiredRVer <- "2.8.0"
    stopifnot(require("utils"))
    if (compareVersion(thisRVer, requiredRVer) < 0) {
        stop(paste("\nYou are currently running R version ", thisRVer, ",\n",
                   "however an R version greater than ", requiredRVer,
                   " is required.\n",
                   "To obtain the latest version of R, please see:\n",
                   "  http://www.r-project.org/\n", sep=""))
    }
    ## update this to reflect most recent supported R version
    if (compareVersion(thisRVer, "3.5.0") >= 0)
        choppedRVer <- "3.5" 
    else
        choppedRVer <- gsub("(\\w+).(\\w+).(\\w+)", "\\1.\\2", thisRVer)
    scriptUrl <- paste(oompaROOT, "installScripts",
                       choppedRVer, "oompaInstall.R", sep="/")
    safeSource <- function() {
        source(scriptUrl, local=TRUE)
        for (objname in ls(all.names=TRUE)) {
            if (exists(objname, envir=.GlobalEnv, inherits=FALSE))
                warning("Redefining ", sQuote(objname))
            .GlobalEnv[[objname]] <- get(objname, inherits=FALSE)
        }
    }
    safeSource()
}

sourceOOMPAinstallScript()

### Install OOMPA packages using CRAN-style repositories.
### ...: arguments passed to install.packages.
getOOMPA <- function(...)
{
    oompainstall(...)
}

oompaLite <- function(pkgs, groupName="lite", ...)
{
    if (missing(pkgs))
        oompainstall(groupName=groupName, ...)
    else
        oompainstall(pkgs=pkgs, groupName=groupName, ...)
}

