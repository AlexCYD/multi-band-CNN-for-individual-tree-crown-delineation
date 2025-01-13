path <- file.path("out", "training", "model")
dirs <- list.dirs(path, recursive = FALSE)

final <- c()
for (numberOfBands in 1:15) {
  selectedDirs <- dirs[grepl(paste0("band", numberOfBands, " "), dirs)]

  out <- c()
  for (selectedDir in selectedDirs) {
    path <- file.path(selectedDir, "trainingHistory.RDS")
    trainingHistory <- readRDS(path)

    # best of epochs
    val_accuracy <- trainingHistory$metrics$val_accuracy
    loss <- trainingHistory$metrics$val_loss
    band <- strsplit(basename(selectedDir), " ")[[1]][2]
    out[band] <- round(val_accuracy[match(min(loss), loss)], 4)
  }

  # best of bands
  out <- out[which.max(out)]
  final[paste0(numberOfBands, "_", names(out))] <- out
}

final
#   1_CHM   2_113  3_ARVI  4_NDVI 5_GNDVI    6_55   7_PRI    8_11     9_B    10_G 
#  0.8554  0.8700  0.8707  0.8701  0.8739  0.8755  0.8739  0.8764  0.8709  0.8667 
#    11_R 12_NDLI  13_EVI  14_NIR 15_SAVI 
#  0.8694  0.8692  0.8627  0.8639  0.8334

# plot
newRes <- 300
f <- newRes / 72
filename <- file.path("out", "analysis", "forward band selection.png")

png(filename, width = 400 * f, height = 200 * f, res = 72 * f, pointsize = 12)
par(mar = c(5, 4, 0, 0) + 0.1, family = "serif")
plot(final, type = "b", xlab = "", ylab = "Validation accuracy", xaxt = "n")
axis(1, at = 1:15, label = names(final), las = 2)
dev.off()
