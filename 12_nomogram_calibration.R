work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(regplot)
library(rms)
out_dir <- file.path(work_dir, "PostLasso", "23.Nomo")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
risk_file <- file.path(work_dir, "PostLasso", "risk.all.txt")
cli_file  <- file.path(work_dir, "clinical-PRAD.txt")
risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table(cli_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- cli[apply(cli, 1, function(x) any(is.na(match('unknow', x)))), , drop = FALSE]
cli$Age <- as.numeric(cli$Age)
cli$PSA <- as.numeric(cli$PSA)
cli$Gleason_Score <- as.numeric(cli$Gleason_Score)
sam_sample <- intersect(row.names(risk), row.names(cli))
risk1 <- risk[sam_sample, , drop = FALSE]
cli <- cli[sam_sample, , drop = FALSE]
rt <- cbind(risk1[, c("futime", "fustat", "risk")], cli)
res_cox <- coxph(Surv(futime, fustat) ~ ., data = rt)
nom1 <- regplot(res_cox,
    plots = c("density", "boxes"),
    clickable = FALSE,
    title = "",
    points = TRUE,
    droplines = TRUE,
    observation = rt[15, ],
    rank = "sd",
    failtime = c(1, 3, 5),
    prfail = FALSE)
nomo_risk <- predict(res_cox, data = rt, type = "risk")
rt <- cbind(risk1, Nomogram = nomo_risk)
out_tab <- rbind(ID = colnames(rt), rt)
write.table(out_tab, file = "nomoRisk.txt", sep = "\t", col.names = FALSE, quote = FALSE)
pdf(file = "calibration.pdf", width = 5, height = 5)
f <- cph(Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 1)
cal <- calibrate(f, cmethod = "KM", method = "boot", u = 1, m = (nrow(rt) / 3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-predicted Survival (%)", ylab = "Observed Survival (%)",
     lwd = 1.5, col = "green", sub = FALSE)
f <- cph(Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 3)
cal <- calibrate(f, cmethod = "KM", method = "boot", u = 3, m = (nrow(rt) / 3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", lwd = 1.5, col = "blue", sub = FALSE, add = TRUE)
f <- cph(Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 5)
cal <- calibrate(f, cmethod = "KM", method = "boot", u = 5, m = (nrow(rt) / 3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", lwd = 1.5, col = "red", sub = FALSE, add = TRUE)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col = c("green", "blue", "red"), lwd = 1.5, bty = 'n')
dev.off()
cat("Module 12 completed\n")
