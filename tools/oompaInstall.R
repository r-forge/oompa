###
### oompaInstall.R
###

## !!! Always change version number when updating this file !!!
oompaVersion <- "3.2"

OOMPAinstallRepos <- function()
{
    ## OOMPA CRAN-style repositories.
    oompa <- "http://silicovore.com/OOMPA"
    setRepositories(FALSE, ind=1:8)
    myRepos <- c(OOMPA = oompa,
                 getOption("repos"))
    myRepos
}

OOMPAinstallPkgGroups <- function(groupName="default")
{
    validGroups <- c(
                     "supercurve",
                     "lite",
                     "tailrank",
                     "genalg",
                     "thresher",
                     "bimodal",
                     "prediction",
                     "default",
                     "arraycube",
                     "all")
    groupName <- match.arg(groupName, validGroups)

    supercurve.group <- c("tclish",
                          "SuperCurve",
                          "SuperCurveGUI",
                          "SlideDesignerGUI",
                          "SuperCurveSampleData")
                          
    lite.group <- c("oompaBase",
                    "oompaData",
                    "PreProcess",
                    "ClassDiscovery",
                    "ClassComparison")

    tailrank.group <- c(lite.group,
                        "TailRank")

    genalg.group <- c(lite.group,
                      "GenAlgo")

    thresher.group <- c(lite.group,
                        "Thresher")

    bimodal.group <- c(lite.group,
                       "BimodalIndex")

    prediction.group <- c(tailrank.group,
                          "Modeler",
                          "CrossValidate",
                          "GenAlgo",
                          "BimodalIndex",
                          "ClassPrediction")

    default.group <- prediction.group
    
    arraycube.group <- c("oompaBase",
                         "ArrayCube",
                         "MINiML")

    switch(EXPR=groupName,
           supercurve = supercurve.group,
           lite = lite.group,
           tailrank = tailrank.group,
           genalg = genalg.group,
           thresher = thresher.group,
           bimodal = bimodal.group,
           prediction = prediction.group,
           arraycube = arraycube.group,
           default = default.group,
           all = {
               OOMPA.url <- OOMPAinstallRepos()[[1]]
               contriburl <- paste(OOMPA.url, "src/contrib",
                                   paste(oompaVersion, "0", sep='.'),
                                   sep="/")
               all.group <- available.packages(contriburl)[, "Package"]
               names(all.group) <- NULL
               all.group
           },
           stop(sprintf("unknown groupName %s", sQuote(groupName))))
}


## OOMPAinstall() version 3.2
## Called by OOMPALite() when R version 3.2 is detected.
##
## Install OOMPA packages using CRAN-style repositories.
## Arguments:
##
##   pkgs: character vector of OOMPA packages to install.
##         The groupName argument will be ignored if pkgs is specified.
##
##   groupName: character matching one of "default", "affy", "graph", "lite",
##              "monograph", "all".
##
##   repos: the user should not try to pass this argument. It's in the argument list
##          of OOMPAinstall() only because we want to raise an error when
##          the user tries to pass it.
##
##   dependencies: passed to install.packages (see install.packages documentation).
##
##   type: passed to install.packages (see install.packages documentation).
##
##   ...: extra arguments passed to install.packages (see install.packages documentation).

oompainstall <- function(pkgs,
                         groupName="default",
                         repos,
                         dependencies=c("Depends", "Imports"),
                         type=getOption("pkgType"),
                         ...)
{
    ## R version 3.2 should have been detected by OOMPALite()
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    cat(paste("Running OOMPAinstall version", oompaVersion, "with R version", thisRVer, "\n"))

    stopifnot(require("utils"))
    cat(paste("Your version of R requires version", oompaVersion, "of OOMPA.\n"))
    
    if (!missing(repos))
        stop("You can't pass a 'repos' argument to the 'OOMPAinstall' function\n")
    repos <- OOMPAinstallRepos()

    if (missing(pkgs)) {
        pkgs <- OOMPAinstallPkgGroups(groupName)
        cat("Will install the following packages:\n")
        print(pkgs)
        cat("Please wait...\n\n")
    }

    install.packages(pkgs=pkgs, repos=repos, dependencies=dependencies, type=type, ...)
}

