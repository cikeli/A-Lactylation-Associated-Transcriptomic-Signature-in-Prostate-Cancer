work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(survival)
library(survminer)

out_dir <- file.path(work_dir, "PostLasso", "PHtest")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

ph_test <- function(input_file, formula_vars, label) {
    rt <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    fml <- as.formula(paste("Surv(futime, fustat) ~", paste(formula_vars, collapse = " + ")))
    fit <- coxph(fml, data = rt)
    zph <- cox.zph(fit)
    res <- cbind(Model = label, as.data.frame(zph$table))
    write.table(res, paste0(label, "_coxzph.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
    pdf(paste0(label, "_coxzph.pdf"), width = 7, height = 6)
    plot(zph)
    dev.off()
    return(res)
}

ph_test(file.path(work_dir, "PostLasso", "lassogeneexp.txt"),
        c("ZNF695", "ATP1A2", "FOXN4"), "TCGA_3gene_model")
ph_test(file.path(work_dir, "PostLasso", "multiCox_input.txt"),
        c("T_stage", "Gleason_Score", "PSA", "riskScore"), "TCGA_multivariable")
cat("Module 22 completed\n")
