work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(reshape2)
library(ggpubr)
out_dir <- file.path(work_dir, "PostLasso", "45.estimateVioplot")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
risk_file <- file.path(work_dir, "risk.all.txt")
tme_file  <- file.path(work_dir, "TMEscores.txt")
risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
risk$risk <- factor(risk$risk, levels = c("low", "high"))
score <- read.table(tme_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
score <- score[, 1:3]
rownames(score) <- gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
score <- score[row.names(risk), , drop = FALSE]
rt <- cbind(risk[, "risk", drop = FALSE], score)
data <- melt(rt, id.vars = c("risk"))
colnames(data) <- c("Risk", "scoreType", "Score")
p <- ggviolin(data, x = "scoreType", y = "Score", fill = "Risk",
    xlab = "",
    ylab = "TME score",
    legend.title = "Risk",
    add = "boxplot", add.params = list(color = "white"),
    palette = c("#0088FF", "#FF5555"), width = 1)
p <- p + rotate_x_text(45)
p1 <- p + stat_compare_means(aes(group = Risk),
    method = "wilcox.test",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
    label = "p.signif")
pdf(file = "vioplot.pdf", width = 6, height = 5)
print(p1)
dev.off()
cat("Module 15 completed\n")
