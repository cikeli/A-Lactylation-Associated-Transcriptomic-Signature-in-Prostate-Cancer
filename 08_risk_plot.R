work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(pheatmap)
out_dir <- file.path(work_dir, "13.riskPlot")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
bio_risk_plot <- function(input_file = NULL, project = NULL) {
    rt <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    rt <- rt[order(rt$riskScore), ]
    risk_class <- rt[, "risk"]
    low_length <- length(risk_class[risk_class == "low"])
    high_length <- length(risk_class[risk_class == "high"])
    low_max <- max(rt$riskScore[risk_class == "low"])
    line <- rt[, "riskScore"]
    line[line > 10] <- 10
    pdf(file = paste0(project, ".riskScore.pdf"), width = 7, height = 4)
    plot(line, type = "p", pch = 20,
         xlab = "Patients (increasing risk score)",
         ylab = "Risk score",
         col = c(rep("blue", low_length), rep("red", high_length)))
    abline(h = low_max, v = low_length, lty = 2)
    legend("topleft", c("High risk", "Low Risk"), bty = "n", pch = 19, col = c("red", "blue"), cex = 1.2)
    dev.off()
    color <- as.vector(rt$fustat)
    color[color == 1] <- "red"
    color[color == 0] <- "blue"
    pdf(file = paste0(project, ".survStat.pdf"), width = 7, height = 4)
    plot(rt$futime, pch = 19,
         xlab = "Patients (increasing risk score)",
         ylab = "Survival time (years)",
         col = color)
    legend("topleft", c("Dead", "Alive"), bty = "n", pch = 19, col = c("red", "blue"), cex = 1.2)
    abline(v = low_length, lty = 2)
    dev.off()
    ann_colors <- list()
    bio_col <- c("blue", "red")
    names(bio_col) <- c("low", "high")
    ann_colors[["Risk"]] <- bio_col
    rt1 <- rt[c(3:(ncol(rt) - 2))]
    rt1 <- t(rt1)
    annotation <- data.frame(Risk = rt[, ncol(rt)])
    rownames(annotation) <- rownames(rt)
    pdf(file = paste0(project, ".heatmap.pdf"), width = 7, height = 4)
    pheatmap(rt1,
        annotation = annotation,
        annotation_colors = ann_colors,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        show_colnames = FALSE,
        scale = "row",
        color = colorRampPalette(c(rep("blue", 3.5), "white", rep("red", 3.5)))(50),
        fontsize_col = 3,
        fontsize = 7,
        fontsize_row = 8)
    dev.off()
}
cat("Module 8 completed\n")
