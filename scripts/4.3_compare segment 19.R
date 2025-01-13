# overlap
n <- 24

# grid search parameters
params <- expand.grid(ws = c(5, 10, 15), hmin = c(0.4, 0.5, 0.6), max_cr = c(30, 50, 70), minCrownsArea = c(0, 1.5, 3, 4.5))

# compare grid search combinations
F1scores <- c()
for (j in seq_len(nrow(params))) {
  path <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j), "imageCrowns_submission.RDS")
  results <- readRDS(path)
  precision <- results$overall$precision
  recall <- results$overall$recall
  F1score <- 2 * precision * recall / (precision + recall)
  F1scores <- c(F1scores, F1score)
}
out <- cbind(params, F1scores)

out[match(max(out$F1scores), out$F1scores), ]
#    ws hmin max_cr minCrownsArea  F1scores
# 76  5  0.5     70             3 0.3759974

summary(out$F1scores)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.1951  0.2669  0.3192  0.3034  0.3469  0.3760

path <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", 76), "imageCrowns_submission.RDS")
results <- readRDS(path)
results$overall$precision
# [1] 0.3527895
results$overall$recall
# [1] 0.4024737
