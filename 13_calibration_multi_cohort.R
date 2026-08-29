work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"
library(survival)
library(ggplot2)
cal_curve <- function(fit, rt, time_point, n_groups = 4) {
    lp <- predict(fit, newdata = rt, type = "lp")
    bf <- basehaz(fit, centered = TRUE)
    H0 <- approx(bf$time, bf$hazard, xout = time_point, rule = 2)$y
    pred_surv <- exp(-H0 * exp(lp))
    groups <- cut(pred_surv,
                  breaks = quantile(pred_surv, probs = seq(0, 1, length.out = n_groups + 1)),
                  include.lowest = TRUE, labels = FALSE)
    result <- data.frame(
        group = 1:n_groups,
        mean_predicted = NA,
        observed_surv = NA,
        n = NA,
        stringsAsFactors = FALSE
    )
    for (g in 1:n_groups) {
        idx <- which(groups == g)
        if (length(idx) < 3) next
        result$n[g] <- length(idx)
        result$mean_predicted[g] <- mean(pred_surv[idx])
        sub_rt <- rt[idx, ]
        km <- survfit(Surv(futime, fustat) ~ 1, data = sub_rt)
        surv_at_t <- approx(km$time, km$surv, xout = time_point, rule = 2)$y
        if (is.na(surv_at_t)) surv_at_t <- min(km$surv)
        result$observed_surv[g] <- surv_at_t
    }
    result <- result[!is.na(result$observed_surv), ]
    return(result)
}
plot_cal <- function(wd, risk_file, label, time_points = c(1, 3, 5)) {
    old_wd <- getwd()
    setwd(wd)
    on.exit(setwd(old_wd))
    risk <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    rt <- risk[, c("futime", "fustat", "riskScore")]
    n <- nrow(rt)
    fit <- coxph(Surv(futime, fustat) ~ riskScore, data = rt)
    cat(sprintf("\n=== %s (n=%d) ===\n", label, n))
    for (t in time_points) {
        cat(sprintf("  %d-year calibration... ", t))
        events_by_t <- sum(rt$fustat == 1 & rt$futime <= t)
        if (events_by_t < 10) {
            cat(sprintf("SKIP (only %d events by year %d)\n", events_by_t, t))
            next
        }
        n_grp <- if (n >= 200) 4 else 3
        cal_df <- cal_curve(fit, rt, time_point = t, n_groups = n_grp)
        cal_df$time_label <- paste0(t, "-Year")
        cat(sprintf("%d groups, %d events\n", nrow(cal_df), events_by_t))
        ax_min <- min(c(cal_df$mean_predicted, cal_df$observed_surv), na.rm = TRUE) - 0.05
        ax_max <- max(c(cal_df$mean_predicted, cal_df$observed_surv), na.rm = TRUE) + 0.05
        ax_min <- max(0, ax_min)
        ax_max <- min(1, ax_max)
        range_val <- ax_max - ax_min
        ax_min <- max(0, ax_min - range_val * 0.1)
        ax_max <- min(1, ax_max + range_val * 0.1)
        p <- ggplot(cal_df, aes(x = mean_predicted, y = observed_surv)) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
            geom_line(color = "#2471A3", linewidth = 1) +
            geom_point(aes(size = n), color = "#2471A3", fill = "#AED6F1", shape = 21, stroke = 1.2) +
            scale_size_continuous(range = c(3, 7), guide = "none") +
            geom_errorbar(aes(
                ymin = pmax(0, observed_surv - 1.96 * sqrt(observed_surv * (1 - observed_surv) / n)),
                ymax = pmin(1, observed_surv + 1.96 * sqrt(observed_surv * (1 - observed_surv) / n))
            ), width = 0.01, color = "#2471A3", alpha = 0.5, linewidth = 0.8) +
            geom_text(aes(label = paste0("n=", n)), vjust = -1.2, size = 3, color = "gray40") +
            labs(x = paste0("Predicted ", t, "-Year Survival Probability"),
                 y = paste0("Observed ", t, "-Year Survival Probability"),
                 title = paste0(label, "  |  ", t, "-Year Calibration")) +
            coord_equal(xlim = c(ax_min, ax_max), ylim = c(ax_min, ax_max)) +
            theme_classic(base_size = 11) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                  panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5))
        pdf(file = paste0("Calibration_", t, "y.pdf"), width = 5.5, height = 5)
        print(p)
        dev.off()
    }
}
base_dir <- file.path(work_dir, "PostLasso", "34.Calibration")
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
cat("Module 13 completed\n")
