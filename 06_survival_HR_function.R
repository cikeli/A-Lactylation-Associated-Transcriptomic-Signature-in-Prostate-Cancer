work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(survminer)
library(writexl)
out_dir <- file.path(work_dir, "15.ROC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
bio_surv <- function(input_file = NULL, surv_pdf = NULL, surv_png = NULL) {
    rt <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)
    rt$risk <- factor(rt$risk, levels = c("low", "high"))
    diff <- survdiff(Surv(futime, fustat) ~ risk, data = rt)
    p_value <- 1 - pchisq(diff$chisq, df = 1)
    if (p_value < 0.001) {
        p_value <- "p<0.001"
    } else {
        p_value <- paste0("p=", sprintf("%.03f", p_value))
    }
    cox <- coxph(Surv(futime, fustat) ~ risk, data = rt)
    cox_sum <- summary(cox)
    HR <- cox_sum$conf.int[1]
    HR_low <- cox_sum$conf.int[3]
    HR_high <- cox_sum$conf.int[4]
    cox_p <- cox_sum$coefficients[5]
    fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
    sur_plot <- ggsurvplot(fit,
        data = rt,
        pval = p_value,
        pval.size = 5,
        pval.coord = c(1, 0.04),
        conf.int = FALSE,
        censor.shape = "|",
        censor.size = 2.5,
        palette = c("blue", "red"),
        linetype = "solid",
        size = 0.6,
        xlab = "Time (years)",
        ylab = "Survival probability",
        legend.title = "",
        legend.labs = c("Low risk", "High risk"),
        legend = c(0.18, 0.22),
        break.time.by = 1,
        font.x = 12,
        font.y = 12,
        font.tickslab = 10,
        font.legend = 10,
        risk.table = TRUE,
        cumevents = TRUE,
        risk.table.col = "strata",
        risk.table.height = 0.25,
        cumevents.height = 0.25,
        risk.table.title = "",
        cumevents.title = "",
        tables.y.text = TRUE,
        ggtheme = theme_classic(),
        surv.median.line = "none")
    pdf(file = surv_pdf, width = 6, height = 7)
    suppressWarnings(suppressMessages(print(sur_plot)))
    dev.off()
    png(file = surv_png, width = 6, height = 7, units = "in", res = 300)
    suppressWarnings(suppressMessages(print(sur_plot)))
    dev.off()
    return(data.frame(
        Dataset = basename(input_file),
        HR = HR,
        HR_low = HR_low,
        HR_high = HR_high,
        Pvalue = cox_p,
        stringsAsFactors = FALSE
    ))
}
cat("Module 6 completed\n")
