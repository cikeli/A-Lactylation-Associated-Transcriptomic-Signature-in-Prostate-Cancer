work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(GSVA)
library(GSEABase)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

exp_file  <- file.path(work_dir, "merge.txt")
risk_file <- file.path(work_dir, "PostLasso", "risk.txt")
out_dir   <- file.path(work_dir, "22.ssGSEA")
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

immune_cell <- read.table(file.path(work_dir, "immuneCell23.txt"), header = FALSE, sep = "\t", check.names = FALSE)
immune_func <- read.table(file.path(work_dir, "immuneFunction13.txt"), header = FALSE, sep = "\t", check.names = FALSE)
cell_sets <- split(immune_cell$V2, immune_cell$V1)
func_sets <- split(immune_func$V2, immune_func$V1)

set.seed(123456)
cell_score <- gsva(data, cell_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
func_score <- gsva(data, func_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
write.table(cell_score, "ssgseaCellScore.txt", sep = "\t", quote = FALSE)
write.table(func_score, "ssgseaFuncScore.txt", sep = "\t", quote = FALSE)

risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
same_sample <- intersect(colnames(cell_score), rownames(risk))
cell_score <- cell_score[, same_sample]
risk <- risk[same_sample, ]

plot_box <- function(score_mat, ylab, outfile) {
    df <- data.frame(t(score_mat), risk = risk[rownames(t(score_mat)), "risk"])
    df <- melt(df, id.vars = "risk")
    colnames(df) <- c("risk", "Term", "Score")
    p <- ggboxplot(df, x = "Term", y = "Score", fill = "risk",
                   palette = c("blue", "red"), xlab = "", ylab = ylab) +
        stat_compare_means(aes(group = risk), method = "wilcox.test",
                           label = "p.signif", label.y = max(df$Score) * 1.02) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    pdf(outfile, width = 9, height = 5)
    print(p)
    dev.off()
}
plot_box(cell_score, "Immune cell infiltration", "ssgseaCell.pdf")
func_score <- func_score[, same_sample]
plot_box(func_score, "Immune function score", "ssgseaFunc.pdf")
cat("Module 19 completed\n")
