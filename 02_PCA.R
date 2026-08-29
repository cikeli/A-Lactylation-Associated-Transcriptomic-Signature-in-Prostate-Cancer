work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(limma)
library(ggplot2)
exp_file    <- file.path(work_dir, "GeneExp.txt")
cluster_file <- file.path(work_dir, "4.Cluster", "Cluster.txt")
out_dir     <- file.path(work_dir, "6.PCA")
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
data <- t(data)
data_pca <- prcomp(data, scale. = TRUE)
pca_predict <- predict(data_pca)
write.table(pca_predict, file = "newTab.xls", quote = FALSE, sep = "\t")
cluster <- read.table(cluster_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cluster <- as.vector(cluster[, 1])
bio_col <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121")
clu_col <- bio_col[1:length(levels(factor(cluster)))]
PCA <- data.frame(PC1 = pca_predict[, 1], PC2 = pca_predict[, 2], cluster = cluster)
PCA_mean <- aggregate(PCA[, 1:2], list(cluster = PCA$cluster), mean)
pdf(file = "PCA.pdf", height = 5, width = 6.5)
ggplot(data = PCA, aes(PC1, PC2)) +
    geom_point(aes(color = cluster)) +
    scale_colour_manual(name = "cluster", values = clu_col) +
    theme_bw() +
    theme(plot.margin = unit(rep(1.5, 4), 'lines')) +
    annotate("text", x = PCA_mean$PC1, y = PCA_mean$PC2, label = PCA_mean$cluster, cex = 7) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
cat("Module 2 completed\n")
