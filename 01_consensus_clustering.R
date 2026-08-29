work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    BiocManager::install("ConsensusClusterPlus")
}
library(ConsensusClusterPlus)
exp_file <- file.path(work_dir, "GeneExp.txt")
out_dir  <- file.path(work_dir, "4.Cluster")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
data <- read.table(exp_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data <- as.matrix(data)
maxK <- 9
results <- ConsensusClusterPlus(
    data,
    maxK = maxK,
    reps = 50,
    pItem = 0.8,
    pFeature = 1,
    title = out_dir,
    clusterAlg = "pam",
    distance = "euclidean",
    seed = 123456,
    plot = "png"
)
cluster_num <- 2
cluster <- results[[cluster_num]][["consensusClass"]]
cluster <- as.data.frame(cluster)
colnames(cluster) <- "cluster"
letter <- c("A", "B", "C", "D", "E", "F", "G")
uniq_clu <- levels(factor(cluster$cluster))
cluster$cluster <- letter[match(cluster$cluster, uniq_clu)]
cluster_out <- rbind(ID = colnames(cluster), cluster)
write.table(cluster_out, file = "Cluster.txt", sep = "\t", quote = FALSE, col.names = FALSE)
cat("Module 1 completed\n")
