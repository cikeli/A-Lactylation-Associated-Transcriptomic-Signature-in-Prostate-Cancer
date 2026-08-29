work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(glmnet)
library(survival)

exp_file <- file.path(work_dir, "PostLasso", "lassoInput.txt")
out_dir  <- file.path(work_dir, "PostLasso", "12.LASSO")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

rt <- read.table(exp_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
x <- as.matrix(rt[, 3:ncol(rt)])
y <- Surv(rt$futime, rt$fustat)

set.seed(123456)
cv_fit <- cv.glmnet(x, y, family = "cox", nfolds = 10, alpha = 1)

pdf("lambda.pdf", width = 6, height = 5)
plot(cv_fit)
abline(v = log(c(cv_fit$lambda.min, cv_fit$lambda.1se)), lty = 2, col = "red")
dev.off()

fit <- glmnet(x, y, family = "cox", alpha = 1)
pdf("coefficient.pdf", width = 6, height = 5)
plot(fit, xvar = "lambda", label = TRUE)
abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
dev.off()

coef_min <- coef(fit, s = cv_fit$lambda.min)
active <- which(coef_min != 0)
gene_coef <- data.frame(Gene = rownames(coef_min)[active],
                        Coefficient = as.numeric(coef_min[active]))
write.table(gene_coef, "lassoCoefficients.txt", sep = "\t", quote = FALSE, row.names = FALSE)

risk_score <- as.numeric(predict(fit, newx = x, s = cv_fit$lambda.min, type = "link"))
risk <- ifelse(risk_score > median(risk_score), "high", "low")
out <- cbind(id = rownames(rt), rt[, c("futime", "fustat")], riskScore = risk_score, risk = risk)
write.table(out, "risk.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Module 18 completed\n")
