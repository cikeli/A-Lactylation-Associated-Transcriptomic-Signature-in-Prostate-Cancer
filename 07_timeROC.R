work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(survminer)
library(timeROC)
out_dir <- file.path(work_dir, "39.ROC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
bio_roc <- function(input_file = NULL, roc_file = NULL) {
    rt <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)
    ROC_rt <- timeROC(
        T = rt$futime,
        delta = rt$fustat,
        marker = rt$riskScore,
        cause = 1,
        weighting = 'aalen',
        times = c(1, 3, 5),
        ROC = TRUE
    )
    pdf(file = roc_file, width = 5, height = 5)
    plot(ROC_rt, time = 1, col = 'green', title = FALSE, lwd = 2)
    plot(ROC_rt, time = 3, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
    plot(ROC_rt, time = 5, col = 'red', add = TRUE, title = FALSE, lwd = 2)
    legend('bottomright',
        c(paste0('AUC at 1 year: ', sprintf("%.03f", ROC_rt$AUC[1])),
          paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
          paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[3]))),
        col = c("green", 'blue', 'red'), lwd = 2, bty = 'n')
    dev.off()
}
cat("Module 7 completed\n")
