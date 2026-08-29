work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(survminer)
out_dir <- file.path(work_dir, "PostLasso", "13.multiCox")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
rt <- read.table(file.path(work_dir, "PostLasso", "lassogeneexp.txt"),
    header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
multi_cox <- coxph(Surv(futime, fustat) ~ ., data = rt)
multi_cox_sum <- summary(multi_cox)
out_tab <- data.frame(
    coef = multi_cox_sum$coefficients[, "coef"],
    HR = multi_cox_sum$conf.int[, "exp(coef)"],
    HR.95L = multi_cox_sum$conf.int[, "lower .95"],
    HR.95H = multi_cox_sum$conf.int[, "upper .95"],
    pvalue = multi_cox_sum$coefficients[, "Pr(>|z|)"]
)
out_tab <- cbind(id = row.names(out_tab), out_tab)
write.table(out_tab, file = "multiCox.xls", sep = "\t", row.names = FALSE, quote = FALSE)
pdf(file = "forest.pdf", width = 8, height = 5)
ggforest(multi_cox,
    main = "Hazard ratio",
    cpositions = c(0.02, 0.22, 0.4),
    fontsize = 0.7,
    refLabel = "reference",
    noDigits = 2)
dev.off()
risk_score <- predict(multi_cox, type = "risk", newdata = rt)
cox_gene <- rownames(multi_cox_sum$coefficients)
cox_gene <- gsub("`", "", cox_gene)
out_col <- c("futime", "fustat", cox_gene)
risk <- as.vector(ifelse(risk_score > median(risk_score), "high", "low"))
write.table(
    cbind(id = rownames(cbind(rt[, out_col], risk_score, risk)),
          cbind(rt[, out_col], risk_score, risk)),
    file = "risk.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
cat("Module 10 completed\n")
