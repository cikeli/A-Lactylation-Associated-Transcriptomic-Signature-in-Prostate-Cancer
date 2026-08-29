work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(ggpubr)
out_dir <- file.path(work_dir, "49.IPS")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
tcia_file  <- file.path(work_dir, "TCIA.txt")
score_file <- file.path(work_dir, "score.group.txt")
ips <- read.table(tcia_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
score <- read.table(score_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
same_sample <- intersect(row.names(ips), row.names(score))
ips <- ips[same_sample, , drop = FALSE]
score <- score[same_sample, "group", drop = FALSE]
data <- cbind(ips, score)
data$group <- ifelse(data$group == "high", "High", "Low")
group <- levels(factor(data$group))
comp <- combn(group, 2)
my_comparisons <- list()
for (i in 1:ncol(comp)) {
    my_comparisons[[i]] <- comp[, i]
}
for (i in colnames(data)[1:(ncol(data) - 1)]) {
    rt <- data[, c(i, "group")]
    colnames(rt) <- c("IPS", "group")
    gg1 <- ggviolin(rt, x = "group", y = "IPS", fill = "group",
        xlab = "Risk Group", ylab = i,
        legend.title = "Risk",
        palette = c("#0066FF", "#FF0000"),
        add = "boxplot", add.params = list(fill = "white")) +
        stat_compare_means(comparisons = my_comparisons)
    pdf(file = paste0(i, ".pdf"), width = 6, height = 5)
    print(gg1)
    dev.off()
}
cat("Module 16 completed\n")
