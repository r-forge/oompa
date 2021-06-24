### sometimes emacs is stupid
plot(2,3)
dev.off()

options(install.packages.compile.from.source = "never")

mypacks <- c("animation",
             "arules",
             "arulesViz",
             "beanplot",
             "beeswarm",
             "caret",
             "clValid",
             "DT",
             "Hmisc",
             "pROC",
             "prodlim",
             "ROCR",
             "rpart.plot",
             "venn",
             "VennDiagram",
             "xlsx",
             "ztable"
             )

for (p in mypacks) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p)
    library(p, character.only=TRUE)
  }
}
