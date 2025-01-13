# sets all random seeds needed to make TensorFlow code reproducible
# this also sets the R random seed set.seed()
tensorflow::set_random_seed(42)

# for GPU also running deterministically
tf$config$experimental$enable_op_determinism()

imageDir <- file.path("out", "training", "image")
maskDir <- file.path("out", "training", "mask")
weightDir <- file.path("out", "training", "weight")
outDir <- file.path("out", "training", "model", "weighted")

model <- unet(c(targetShape, 8), weighted = TRUE)
optimizer <- keras::optimizer_sgd(learning_rate = 0.001, momentum = 0.99)
loss <- myLoss_mean()
metrics <- "accuracy"

epochs <- 50
callbacks <- list(
  keras::callback_early_stopping(
    monitor = "val_loss",
    min_delta = 1e-4,
    patience = 5,
    verbose = 1,
    restore_best_weights = TRUE
  ),
  keras::callback_model_checkpoint(
    file.path(outDir, "checkpoint.ckpt"),
    verbose = 1,
    save_weights_only = TRUE,
    save_best_only = TRUE
  )
)

sampling <- sample(seq_along(list.files(maskDir)), 100)
selectedBands <- readRDS(file.path("out", "analysis", "selectedBands.RDS"))

trainingHistory <- myTrain(
  imageDir, maskDir, weightDir, outDir,
  model, optimizer, loss, metrics,
  batch_size, epochs, callbacks,
  class_weight = NULL, sampling, augment = TRUE, weighted = TRUE,
  useBands = selectedBands[1:8]
)
