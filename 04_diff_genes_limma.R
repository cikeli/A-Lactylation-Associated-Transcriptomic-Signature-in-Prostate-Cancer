work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(limma)
library(VennDiagram)
exp_file    <- file.path(work_dir, "merge.txt")
clu_file    <- file.path(work_dir, "4.Cluster", "Cluster.txt")
out_dir     <- file.path(work_dir, "7.clusterDiff")
adj_p_filter <- 0.05
logfc_filter <- 1
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
cluster <- read.table(clu_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
same_sample <- intersect(colnames(data), row.names(cluster))
data <- data[, same_sample]
cluster <- cluster[same_sample, ]
gene_list <- list()
type <- as.vector(cluster)
design <- model.matrix(~ 0 + factor(type))
colnames(design) <- levels(factor(type))
comp <- combn(levels(factor(type)), 2)
all_diff_genes <- c()
for (i in 1:ncol(comp)) {
    fit <- lmFit(data, design)
    contrast <- paste0(comp[2, i], "-", comp[1, i])
    cont.matrix <- makeContrasts(contrasts = contrast, levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    all_diff <- topTable(fit2, adjust = 'fdr', number = 200000)
    all_diff_out <- rbind(id = colnames(all_diff), all_diff)
    write.table(all_diff_out, file = paste0(contrast, ".all.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    diff_sig <- all_diff[with(all_diff, (abs(logFC) > logfc_filter & P.Value < adj_p_filter)), ]
    diff_sig_out <- rbind(id = colnames(diff_sig), diff_sig)
    write.table(diff_sig_out, file = paste0(contrast, ".diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    gene_list[[contrast]] <- row.names(diff_sig)
}
venn.plot <- venn.diagram(gene_list, filename = NULL, fill = rainbow(length(gene_list)))
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()
inter_genes <- Reduce(intersect, gene_list)
write.table(file = "interGene.txt", inter_genes, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
inter_gene_exp <- data[inter_genes, ]
inter_gene_exp <- rbind(id = colnames(inter_gene_exp), inter_gene_exp)
write.table(inter_gene_exp, file = "interGeneExp.txt", sep = "\t", quote = FALSE, col.names = FALSE)
cat("Module 4 completed\n")
