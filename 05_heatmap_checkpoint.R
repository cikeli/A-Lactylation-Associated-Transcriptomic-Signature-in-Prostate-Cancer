work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(limma)
library(pheatmap)
exp_file   <- file.path(work_dir, "checkpoint.txt")
risk_file  <- file.path(work_dir, "geneCluster.txt")
out_dir    <- file.path(work_dir, "heatmap")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
rt <- read.table(exp_file, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
mat <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
mat <- avereps(mat)
mat <- mat[rowMeans(mat) > 0, ]
data <- t(avereps(mat))
risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
low_sample <- row.names(risk[risk$geneCluster == "A", ])
high_sample <- row.names(risk[risk$geneCluster == "B", ])
low_data <- data[, low_sample]
high_data <- data[, high_sample]
data <- cbind(low_data, high_data)
con_num <- ncol(low_data)
treat_num <- ncol(high_data)
sample_type <- c(rep(1, con_num), rep(2, treat_num))
sig_vec <- c()
for (i in row.names(data)) {
    test <- wilcox.test(data[i, ] ~ sample_type)
    pvalue <- test$p.value
    sig <- ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
    sig_vec <- c(sig_vec, paste0(i, sig))
}
row.names(data) <- sig_vec
type <- c(rep("A", con_num), rep("B", treat_num))
type <- factor(type, levels = c("A", "B"))
names(type) <- colnames(data)
type <- as.data.frame(type)
bio_col <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121")
ann_colors <- list()
clu_col <- bio_col[1:length(levels(factor(type$type)))]
names(clu_col) <- levels(factor(type$type))
ann_colors[["Type"]] <- clu_col
pdf("heatmap-checkpoint.pdf", width = 6, height = 5)
pheatmap(data,
    annotation = type,
    color = colorRampPalette(c(rep("blue", 5), "white", rep("red", 5)))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_colors = ann_colors,
    scale = "row",
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize = 7,
    fontsize_row = 7,
    fontsize_col = 7)
dev.off()
cat("Module 5 completed\n")
