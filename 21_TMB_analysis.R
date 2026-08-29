work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(maftools)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

maf_file <- file.path(work_dir, "TCGA.PRAD.mutect.maf")
risk_file <- file.path(work_dir, "PostLasso", "risk.txt")
cli_file  <- file.path(work_dir, "time.txt")
out_dir   <- file.path(work_dir, "30.TMB")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

maf <- read.maf(maf_file)
tmb <- tmb(maf, captureSize = 38, logScale = FALSE)
tmb$Tumor_Sample_Barcode <- substr(tmb$Tumor_Sample_Barcode, 1, 12)

risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli  <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

rt <- merge(data.frame(id = tmb$Tumor_Sample_Barcode, TMB = tmb$total_perMB),
            data.frame(id = rownames(risk), riskScore = risk$riskScore, risk = risk$risk), by = "id")
rt <- merge(rt, data.frame(id = rownames(cli), futime = cli$futime, fustat = cli$fustat), by = "id")
write.table(rt, "TMB_riskScore.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("cor_TMB_LRS.pdf", width = 5, height = 5)
print(ggscatter(rt, x = "riskScore", y = "TMB", add = "reg.line",
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Risk score", ylab = "TMB (mutations/Mb)"))
dev.off()

pdf("box_TMB.pdf", width = 5, height = 5)
print(ggboxplot(rt, x = "risk", y = "TMB", fill = "risk", palette = c("blue", "red"),
                xlab = "", ylab = "TMB (mutations/Mb)") +
      stat_compare_means(method = "wilcox.test"))
dev.off()

rt$TMB_group <- ifelse(rt$TMB > median(rt$TMB), "high", "low")
run_km <- function(data, group_col, legend_title, outfile) {
    fit <- survfit(Surv(futime, fustat) ~ data[[group_col]], data = data)
    p <- ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE,
                    palette = c("blue", "red"), xlab = "Time (years)",
                    legend.title = legend_title, legend.labs = c("High", "Low"))
    pdf(outfile, width = 6, height = 6.5)
    print(p)
    dev.off()
}
run_km(rt, "TMB_group", "TMB", "KM_TMB.pdf")

rt$combine <- paste(rt$TMB_group, rt$risk, sep = "_")
fit <- survfit(Surv(futime, fustat) ~ combine, data = rt)
pdf("KM_TMB_LRS.pdf", width = 6.5, height = 6.5)
print(ggsurvplot(fit, data = rt, pval = TRUE, risk.table = TRUE,
                 palette = c("#0066FF", "#FF9900", "#6E568C", "#D20A13"),
                 xlab = "Time (years)", legend.title = "TMB + LRS"))
dev.off()

cox_int <- coxph(Surv(futime, fustat) ~ riskScore * TMB, data = rt)
print(summary(cox_int))
sink("interaction_LRS_TMB.txt")
print(summary(cox_int))
sink()
cat("Module 21 completed\n")
