selectedDir <- file.path("out", "training", "model", "weighted")
path <- file.path(selectedDir, "trainingHistory.RDS")
trainingHistory <- readRDS(path)

val_accuracy <- trainingHistory$metrics$val_accuracy
loss <- trainingHistory$metrics$val_loss
round(val_accuracy[match(min(loss), loss)], 4)
# [1] 0.8539
