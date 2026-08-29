work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(ggalluvial)
library(ggplot2)
library(dplyr)
out_dir <- file.path(work_dir, "PostLasso", "39.ggalluvial")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
clu_file     <- file.path(work_dir, "4.Cluster", "Cluster.txt")
gene_clu_file <- file.path(work_dir, "PostLasso", "geneCluster.txt")
score_file   <- file.path(work_dir, "score.group.txt")
cli_file     <- file.path(work_dir, "time.txt")
trait_col    <- "fustat"
clu <- read.table(clu_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
gene_clu <- read.table(gene_clu_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
score <- read.table(score_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
two_cluster <- cbind(clu, gene_clu)
rownames(two_cluster) <- gsub("(.*?)\\_(.*?)", "\\2", rownames(two_cluster))
same_sample <- intersect(row.names(two_cluster), row.names(score))
score_clu <- cbind(score[same_sample, , drop = FALSE], two_cluster[same_sample, , drop = FALSE])
same_sample <- intersect(row.names(score_clu), row.names(cli))
rt <- cbind(score_clu[same_sample, ], cli[same_sample, ])
rt <- rt[, c("cluster", "geneCluster", "group", trait_col)]
colnames(rt) <- c("cluster", "geneCluster", "score", trait_col)
cor_lodes <- to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")
pdf(file = "ggalluvial.pdf", width = 6, height = 6)
mycol <- rep(c("#0066FF", "#FF9900", "#FF0000", "#029149", "#6E568C", "#E0367A", "#D8D155", "#223D6C", "#D20A13", "#431A3D"), 2)
ggplot(cor_lodes, aes(x = x, stratum = stratum, alluvium = Cohort, fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0)) +
    geom_flow(width = 2 / 10, aes.flow = "forward") +
    geom_stratum(alpha = 0.9, width = 2 / 10) +
    scale_fill_manual(values = mycol) +
    geom_text(stat = "stratum", size = 3, color = "black") +
    xlab("") + ylab("") + theme_bw() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank()) +
    ggtitle("") + guides(fill = FALSE)
dev.off()
cat("Module 14 completed\n")
