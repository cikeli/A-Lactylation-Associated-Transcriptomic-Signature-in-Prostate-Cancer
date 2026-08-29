work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(survminer)
cluster_file <- file.path(work_dir, "4.Cluster", "Cluster.txt")
cli_file     <- file.path(work_dir, "time.txt")
out_dir      <- file.path(work_dir, "5.ClusterSur")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
cluster <- read.table(cluster_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rownames(cluster) <- gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
colnames(cli) <- c("futime", "fustat")
same_sample <- intersect(row.names(cluster), row.names(cli))
rt <- cbind(cli[same_sample, , drop = FALSE], cluster[same_sample, , drop = FALSE])
length <- length(levels(factor(rt$cluster)))
diff <- survdiff(Surv(futime, fustat) ~ cluster, data = rt)
p_value <- 1 - pchisq(diff$chisq, df = length - 1)
if (p_value < 0.001) {
    p_value <- "p<0.001"
} else {
    p_value <- paste0("p=", sprintf("%.03f", p_value))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)
bio_col <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121")
bio_col <- bio_col[1:length]
sur_plot <- ggsurvplot(fit,
    data = rt,
    conf.int = FALSE,
    pval = p_value,
    pval.size = 6,
    legend.title = "lacCluster",
    legend.labs = levels(factor(rt[, "cluster"])),
    legend = c(0.8, 0.8),
    font.legend = 10,
    xlab = "Time (years)",
    break.time.by = 1,
    palette = bio_col,
    surv.median.line = "hv",
    risk.table = TRUE,
    cumevents = FALSE,
    risk.table.height = 0.25)
pdf(file = "survival.pdf", onefile = FALSE, width = 7, height = 5.5)
print(sur_plot)
dev.off()
cat("Module 3 completed\n")
