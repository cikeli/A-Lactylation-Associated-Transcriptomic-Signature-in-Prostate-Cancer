work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(plyr)
library(ggplot2)
library(ggpubr)
out_dir <- file.path(work_dir, "40.scoreCli")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
score_file <- file.path(work_dir, "ICIscore.group.txt")
cli_file   <- file.path(work_dir, "clinical-PRAD.txt")
trait_col  <- "T_stage"
score <- read.table(score_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
same_sample <- intersect(row.names(score), row.names(cli))
rt <- cbind(score[same_sample, ], cli[same_sample, ])
rt1 <- rt[, c(trait_col, "group")]
colnames(rt1) <- c("trait", "group")
df <- as.data.frame(table(rt1))
df <- ddply(df, .(group), transform, percent = Freq / sum(Freq) * 100)
df <- ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label <- paste0(sprintf("%.0f", df$percent), "%")
bio_col <- c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121")
bio_col <- bio_col[1:length(unique(rt[, trait_col]))]
p <- ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
    geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
    scale_fill_manual(values = bio_col) +
    xlab("Risk Group") + ylab("Percent") +
    guides(fill = guide_legend(title = trait_col)) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
    theme_bw()
pdf(file = "barplot.pdf", width = 4, height = 5)
print(p)
dev.off()
rt2 <- rt[, c(trait_col, "score")]
colnames(rt2) <- c("trait", "score")
type <- levels(factor(rt2[, "trait"]))
comp <- combn(type, 2)
my_comparisons <- list()
for (i in 1:ncol(comp)) {
    my_comparisons[[i]] <- comp[, i]
}
boxplot <- ggboxplot(rt2, x = "trait", y = "score", fill = "trait",
    xlab = trait_col,
    ylab = "Risk Score",
    legend.title = trait_col,
    palette = bio_col
) + stat_compare_means(comparisons = my_comparisons)
pdf(file = "boxplot.pdf", width = 4, height = 4.5)
print(boxplot)
dev.off()
cat("Module 9 completed\n")
