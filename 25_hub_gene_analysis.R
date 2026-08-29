work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(limma)
library(survival)
library(survminer)
library(pROC)
library(ggplot2)
library(ggpubr)

exp_file  <- file.path(work_dir, "merge.txt")
cli_file  <- file.path(work_dir, "time.txt")
type_file <- file.path(work_dir, "sampleType.txt")
hub_genes <- c("ZNF695", "ATP1A2", "FOXN4")
out_dir   <- file.path(work_dir, "18.hubGene")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

rt <- read.table(exp_file, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]

stype <- read.table(type_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

for (g in hub_genes) {
    expr <- data[g, ]
    df <- data.frame(id = colnames(data), expr = as.numeric(expr),
                     type = stype[colnames(data), 1])
    pdf(paste0(g, "_expression.pdf"), width = 5, height = 5)
    print(ggboxplot(df, x = "type", y = "expr", fill = "type",
                    palette = c("blue", "red"), xlab = "", ylab = paste(g, "expression")) +
          stat_compare_means(method = "wilcox.test"))
    dev.off()

    same_sample <- intersect(names(expr), rownames(cli))
    sf <- data.frame(futime = cli[same_sample, "futime"], fustat = cli[same_sample, "fustat"],
                     expr = as.numeric(expr[same_sample]))
    sf$group <- ifelse(sf$expr > median(sf$expr), "high", "low")
    fit <- survfit(Surv(futime, fustat) ~ group, data = sf)
    pdf(paste0(g, "_survival.pdf"), width = 5.5, height = 6)
    print(ggsurvplot(fit, data = sf, pval = TRUE, risk.table = TRUE,
                     palette = c("blue", "red"), xlab = "Time (years)",
                     legend.title = g, legend.labs = c("High", "Low")))
    dev.off()

    roc_obj <- roc(sf$fustat, sf$expr, quiet = TRUE, ci = TRUE)
    pdf(paste0(g, "_ROC.pdf"), width = 5, height = 5)
    plot(roc_obj, col = "red", lwd = 2, main = paste(g, "ROC"))
    text(0.4, 0.2, paste0("AUC = ", sprintf("%.3f", as.numeric(auc(roc_obj)))))
    dev.off()
}
cat("Module 25 completed\n")
