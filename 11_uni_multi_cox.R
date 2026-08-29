work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(survminer)
out_dir <- file.path(work_dir, "PostLasso", "20.indep")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
bio_forest <- function(cox_file = NULL, forest_file = NULL, forest_col = NULL) {
    rt <- read.table(cox_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    gene <- rownames(rt)
    hr <- sprintf("%.3f", rt$"HR")
    hr_low <- sprintf("%.3f", rt$"HR.95L")
    hr_high <- sprintf("%.3f", rt$"HR.95H")
    hazard_ratio <- paste0(hr, "(", hr_low, "-", hr_high, ")")
    p_val <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    pdf(file = forest_file, width = 6.6, height = 4.5)
    n <- nrow(rt)
    n_row <- n + 1
    ylim <- c(1, n_row)
    layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
    xlim <- c(0.1, 3)
    par(mar = c(4, 2.5, 2, 1))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
    text.cex <- 0.8
    text(0, n:1, gene, adj = 0, cex = text.cex)
    text(1.5 - 0.5 * 0.2, n:1, p_val, adj = 1, cex = text.cex)
    text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)
    text(3.1, n:1, hazard_ratio, adj = 1, cex = text.cex)
    text(3.1, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
    par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
    xlim <- c(0.5, 3)
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
    arrows(as.numeric(hr_low), n:1, as.numeric(hr_high), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
    abline(v = 1, col = "black", lty = 2, lwd = 2)
    boxcolor <- ifelse(as.numeric(hr) > 1, forest_col, forest_col)
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
    axis(1)
    dev.off()
}
indep <- function(risk_file = NULL, cli_file = NULL, project = NULL) {
    risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    cli <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    same_sample <- intersect(row.names(cli), row.names(risk))
    risk <- risk[same_sample, ]
    cli <- cli[same_sample, ]
    rt <- cbind(futime = risk[, 1], fustat = risk[, 2], cli, riskScore = risk[, (ncol(risk) - 1)])
    uni_cox_file <- paste0(project, ".uniCox.txt")
    uni_cox_pdf  <- paste0(project, ".uniCox.pdf")
    uni_tab <- data.frame()
    for (i in colnames(rt[, 3:ncol(rt)])) {
        cox <- coxph(Surv(futime, fustat) ~ rt[, i], data = rt)
        cox_summary <- summary(cox)
        uni_tab <- rbind(uni_tab,
            cbind(id = i,
                  HR = cox_summary$conf.int[, "exp(coef)"],
                  HR.95L = cox_summary$conf.int[, "lower .95"],
                  HR.95H = cox_summary$conf.int[, "upper .95"],
                  pvalue = cox_summary$coefficients[, "Pr(>|z|)"]))
    }
    write.table(uni_tab, file = uni_cox_file, sep = "\t", row.names = FALSE, quote = FALSE)
    bio_forest(cox_file = uni_cox_file, forest_file = uni_cox_pdf, forest_col = "green")
    multi_cox_file <- paste0(project, ".multiCox.txt")
    multi_cox_pdf  <- paste0(project, ".multiCox.pdf")
    uni_tab <- uni_tab[as.numeric(uni_tab[, "pvalue"]) < 1, ]
    rt1 <- rt[, c("futime", "fustat", as.vector(uni_tab[, "id"]))]
    multi_cox <- coxph(Surv(futime, fustat) ~ ., data = rt1)
    multi_cox_sum <- summary(multi_cox)
    multi_tab <- data.frame(
        HR = multi_cox_sum$conf.int[, "exp(coef)"],
        HR.95L = multi_cox_sum$conf.int[, "lower .95"],
        HR.95H = multi_cox_sum$conf.int[, "upper .95"],
        pvalue = multi_cox_sum$coefficients[, "Pr(>|z|)"])
    multi_tab <- cbind(id = row.names(multi_tab), multi_tab)
    write.table(multi_tab, file = multi_cox_file, sep = "\t", row.names = FALSE, quote = FALSE)
    bio_forest(cox_file = multi_cox_file, forest_file = multi_cox_pdf, forest_col = "red")
}
cat("Module 11 completed\n")
