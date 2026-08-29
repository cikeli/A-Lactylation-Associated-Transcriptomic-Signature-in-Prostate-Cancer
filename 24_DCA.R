work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(survival)
library(dcurves)
library(ggplot2)

out_dir <- file.path(work_dir, "PostLasso", "DCA")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

run_dca <- function(input_file, time_point, label) {
    rt <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    d <- dca(Surv(futime, fustat) ~ riskScore, data = rt,
             time = time_point, thresholds = seq(0, 0.5, by = 0.01))
    p <- plot(d, smooth = TRUE) +
        labs(title = paste0(label, " - ", time_point, "-year DCA")) +
        theme_classic()
    pdf(paste0("DCA_", label, "_", time_point, "y.pdf"), width = 6, height = 5)
    print(p)
    dev.off()
}

run_dca(file.path(work_dir, "PostLasso", "risk.txt"), 1, "TCGA")
run_dca(file.path(work_dir, "PostLasso", "risk.txt"), 3, "TCGA")
cat("Module 24 completed\n")
