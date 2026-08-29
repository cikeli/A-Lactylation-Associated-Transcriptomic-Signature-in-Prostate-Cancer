work_dir <- "C:/Users/cikeli/Desktop/YSL-LACTATE"

library(pROC)

out_dir <- file.path(work_dir, "PostLasso", "DeLong")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

time_point <- 1

rt <- read.table(file.path(work_dir, "PostLasso", "rocCompare_input.txt"),
                 header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rt <- rt[rt$futime > time_point | rt$fustat == 1, ]
rt$status <- ifelse(rt$fustat == 1 & rt$futime <= time_point, 1, 0)

markers <- setdiff(colnames(rt), c("futime", "fustat", "status"))
roc_list <- list()
for (m in markers) {
    roc_list[[m]] <- roc(rt$status, rt[[m]], quiet = TRUE, ci = TRUE)
}

aucs <- sapply(roc_list, function(r) as.numeric(auc(r)))
print(round(aucs, 3))

pairs <- combn(markers, 2, simplify = FALSE)
res <- data.frame()
for (pr in pairs) {
    tt <- roc.test(roc_list[[pr[1]]], roc_list[[pr[2]]], method = "delong", paired = TRUE)
    res <- rbind(res, data.frame(Marker1 = pr[1], Marker2 = pr[2],
                                 AUC1 = round(aucs[pr[1]], 3), AUC2 = round(aucs[pr[2]], 3),
                                 Difference = round(aucs[pr[1]] - aucs[pr[2]], 3),
                                 P_value = signif(tt$p.value, 3)))
}
write.table(res, "DeLong_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
cat("Module 23 completed\n")
