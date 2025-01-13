dirList <- c(
  file.path("out", "training", "model", "level1"),
  file.path("out", "training", "model", "level2"),
  file.path("out", "training", "model", "level3"),
  file.path("out", "training", "model", "band8 11"),
  file.path("out", "training", "model", "level5")
)

final <- c()
for (i in dirList) {
  path <- file.path(i, "trainingHistory.RDS")
  trainingHistory <- readRDS(path)

  val_accuracy <- trainingHistory$metrics$val_accuracy
  loss <- trainingHistory$metrics$val_loss
  final[basename(i)] <- round(val_accuracy[match(min(loss), loss)], 4)
}

final
#   level1   level2   level3 band8 11   level5 
#   0.8629   0.8668   0.8688   0.8764   0.8684

# plot
newRes <- 300
f <- newRes / 72

png(file.path("out", "analysis", "compare level.png"), width = 400 * f, height = 200 * f, res = 72 * f, pointsize = 12)
par(mar = c(4, 4, 0, 0) + 0.1, family = "serif")
plot(final, type = "b", xlab = "Number of levels", ylab = "Validation accuracy", yaxt = "n")
axis(2, at = c(0.866, 0.870, 0.874))
dev.off()
