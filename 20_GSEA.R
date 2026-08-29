work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)

exp_file <- file.path(work_dir, "merge.txt")
out_dir  <- file.path(work_dir, "25.GSEA")
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

target_gene <- "ZNF695"
group <- ifelse(data[target_gene, ] > median(data[target_gene, ]), "high", "low")
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- c("high", "low")
fit <- lmFit(data, design)
cont <- makeContrasts("high-low", levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
deg <- topTable(fit2, number = Inf)

gene_list <- deg$t
names(gene_list) <- rownames(deg)
gene_list <- sort(gene_list, decreasing = TRUE)

id_map <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- gene_list[id_map$SYMBOL]
names(gene_list) <- id_map$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

set.seed(123456)
gsea <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "ALL",
              pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.table(as.data.frame(gsea), "gseaGO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

kegg <- gseKEGG(geneList = gene_list, organism = "hsa",
                pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.table(as.data.frame(kegg), "gseaKEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

hallmark <- read.gmt(file.path(work_dir, "h.all.v2023.2.Hs.symbols.gmt"))
gene_list_sym <- deg$t
names(gene_list_sym) <- rownames(deg)
gene_list_sym <- sort(gene_list_sym, decreasing = TRUE)
gsea_h <- GSEA(gene_list_sym, TERM2GENE = hallmark, pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.table(as.data.frame(gsea_h), "gseaHallmark.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("gsea_top5.pdf", width = 8, height = 6)
print(gseaplot2(gsea_h, 1:min(5, nrow(gsea_h)), pvalue_table = TRUE))
dev.off()
cat("Module 20 completed\n")
