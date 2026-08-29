work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(survival)

out_dir <- file.path(work_dir, "PostLasso", "PHtest")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

cohorts <- list(
    Taylor = file.path(work_dir, "taylor", "risk.txt"),
    ICGC   = file.path(work_dir, "ICGC", "risk.txt")
)
for (nm in names(cohorts)) {
    rt <- read.table(cohorts[[nm]], header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    fit <- coxph(Surv(futime, fustat) ~ riskScore, data = rt)
    zph <- cox.zph(fit)
    print(nm)
    print(zph$table)
    write.table(as.data.frame(zph$table), paste0(nm, "_coxzph.txt"), sep = "\t", quote = FALSE)
    pdf(paste0(nm, "_coxzph.pdf"), width = 6, height = 5)
    plot(zph)
    dev.off()
}
cat("Module 26 completed\n")
