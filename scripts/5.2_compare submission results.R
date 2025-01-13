# read submission results
submissions <- list.files(file.path("out", "submission"), full.name = TRUE)
for (submission in submissions) {
  f <- tools::file_path_sans_ext(basename(submission))
  path <- file.path("out", "submission", paste0(f, ".RDS"))
  results <- readRDS(path)
  assign(substr(f, 13, nchar(f)), results)
}

# read my results
path <- file.path("out", "testing", "segment", "band8 11 overlap24 segment76", "imageCrowns_submissionCommon.RDS")
This_study <- readRDS(path)

# overall
overall <- rbind(This_study$overall, Dalponte2016$overall, Li2012$overall, Silva2016$overall, Weinstein_unpublished$overall, Weinstein2019$overall)
submission <- c("This_study", "Dalponte2016", "Li2012", "Silva2016", "Weinstein_unpublished", "Weinstein2019")
overall <- cbind(submission, overall)

overall
#              submission  precision    recall
# 1            This_study 0.40988889 0.5775873
# 2          Dalponte2016 0.31566667 0.5972857
# 3                Li2012 0.07788889 0.2659365
# 4             Silva2016 0.32352381 0.5822222
# 5 Weinstein_unpublished 0.71807937 0.8287778
# 6         Weinstein2019 0.62993651 0.7976508

# plot precision vs recall
newRes <- 300
f <- newRes / 72
x <- overall$precision
y <- overall$recall

png(file.path("out", "analysis", "compare submissions.png"), width = 400 * f, height = 400 * f, res = 72 * f, pointsize = 12)
par(mar = c(4, 4, 0, 0) + 0.1, family = "serif")
plot(x, y, xlim = c(0, 1), ylim = c(0, 1), xlab = "Precision", ylab = "Recall")
text(x[1], y[1], submission[1], pos = 4)
text(x[2], y[2] + 0.01, submission[2], pos = 3)
text(x[3], y[3], submission[3], pos = 4)
text(x[4], y[4] - 0.01, submission[4], pos = 1)
text(x[5], y[5] + 0.005, submission[5], pos = 2)
text(x[6], y[6] - 0.005, submission[6], pos = 2)
dev.off()
