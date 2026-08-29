work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
out_dir <- file.path(work_dir, "PostLasso", "49.pRRophetic")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
set.seed(12345)
p_filter <- 0.001
exp_file <- file.path(work_dir, "TCGA.TPM.txt")
risk_file <- file.path(work_dir, "risk.all.txt")
data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
all_drugs <- unique(drugData2016$Drug.name)
rt <- read.table(exp_file, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0.5, ]
group <- sapply(strsplit(colnames(data), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
data <- data[, group == 0]
data <- t(data)
rownames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data <- avereps(data)
data <- t(data)
risk_rt <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
risk_rt$riskScore[risk_rt$riskScore > quantile(risk_rt$riskScore, 0.99)] <- quantile(risk_rt$riskScore, 0.99)
for (drug in all_drugs) {
    possible_error <- tryCatch({
        sensitivity <- pRRopheticPredict(data, drug, selection = 1, dataset = "cgp2016")
    }, error = function(e) e)
    if (inherits(possible_error, "error")) next
    sensitivity <- sensitivity[sensitivity != "NaN"]
    sensitivity[sensitivity > quantile(sensitivity, 0.99)] <- quantile(sensitivity, 0.99)
    same_sample <- intersect(row.names(risk_rt), names(sensitivity))
    risk <- risk_rt[same_sample, c("riskScore", "risk"), drop = FALSE]
    sensitivity <- sensitivity[same_sample]
    rt <- cbind(risk, sensitivity)
    rt$risk <- factor(rt$risk, levels = c("low", "high"))
    type <- levels(factor(rt[, "risk"]))
    comp <- combn(type, 2)
    my_comparisons <- list()
    for (i in 1:ncol(comp)) {
        my_comparisons[[i]] <- comp[, i]
    }
    test <- wilcox.test(sensitivity ~ risk, data = rt)
    diff_pvalue <- test$p.value
    x <- as.numeric(rt[, "riskScore"])
    y <- as.numeric(rt[, "sensitivity"])
    cor_t <- cor.test(x, y, method = "spearman")
    cor_pvalue <- cor_t$p.value
    if ((diff_pvalue < p_filter) & (cor_pvalue < p_filter)) {
        boxplot <- ggboxplot(rt, x = "risk", y = "sensitivity", fill = "risk",
            xlab = "Risk",
            ylab = paste0(drug, " sensitivity (IC50)"),
            legend.title = "Risk",
            palette = c("#0066FF", "#FF0000")
        ) + stat_compare_means(comparisons = my_comparisons)
        pdf(file = paste0("drugSensitivity.", drug, ".pdf"), width = 5, height = 4.5)
        print(boxplot)
        dev.off()
        df1 <- as.data.frame(cbind(x, y))
        p1 <- ggplot(df1, aes(x, y)) +
            xlab("Risk score") + ylab(paste0(drug, " sensitivity (IC50)")) +
            geom_point() + geom_smooth(method = "lm", formula = y ~ x) + theme_bw() +
            stat_cor(method = 'spearman', aes(x = x, y = y))
        pdf(file = paste0("Cor.", drug, ".pdf"), width = 5, height = 4.6)
        print(p1)
        dev.off()
    }
}
cat("Module 17 completed\n")
cat("Done\n")
