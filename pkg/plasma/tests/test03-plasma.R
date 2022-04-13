library(plasma)
data("TCGA-ESCA")
assemble <- assemble[-1] # remove clinical data until we convert to numeric
MO <- prepareMultiOmics(assemble, Outcome)
## test complete cox models
bigfit <- fitCoxModels(MO, "Days", "vital_status", "dead")
class(bigfit)
summary(bigfit)

extension <- extendCoxModels(MO, bigfit)
class(extension)

bonkers <- data.frame(MO@outcome[, c(2,5)], extension)

model <- coxph(Surv(Days, vital_status == "dead") ~ ., bonkers)
model1 <- step(model, trace = 0)
model1

bonkers$Risk <- predict(model1)
lh <- c("Low", "High")
bonkers$Split <- factor(lh[1 + 1*(bonkers$Risk >- median(bonkers$Risk))],
                        levels = lh)

smodel <- coxph(Surv(Days, vital_status == "dead") ~ Split, bonkers)
plot(survfit(Surv(Days, vital_status == "dead") ~ Split, bonkers),
     col = c("blue", "red"), lwd=2, mark.time = TRUE)
legend("topright", paste(lh, "Risk"), lwd=2, col = c("blue", "red"))

