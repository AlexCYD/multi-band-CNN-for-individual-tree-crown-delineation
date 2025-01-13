bands <- c("R", "G", "B", "CHM", "11", "55", "113", "NIR", "NDVI", "EVI", "ARVI", "PRI", "NDLI", "SAVI", "GNDVI")
selected <- c()
bandsToGo <- 1

for (numberOfBands in bandsToGo:15) {
  out <- c()

  for (band in bands) {
    # sets all random seeds needed to make TensorFlow code reproducible
    # this also sets the R random seed set.seed()
    tensorflow::set_random_seed(42)

    # for GPU also running deterministically
    tf$config$experimental$enable_op_determinism()

    imageDir <- file.path("out", "training", "image")
    maskDir <- file.path("out", "training", "mask")
    outDir <- file.path("out", "training", "model", paste0("band", numberOfBands, " ", band))

    model <- unet(c(targetShape, numberOfBands))
    optimizer <- keras::optimizer_sgd(learning_rate = 0.001, momentum = 0.99)
    loss <- "binary_crossentropy"
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

    trainingHistory <- myTrain(
      imageDir, maskDir, weightDir, outDir,
      model, optimizer, loss, metrics,
      batch_size, epochs, callbacks,
      class_weight = NULL, sampling, augment = TRUE, weighted = FALSE,
      useBands = c(selected, band)
    )

    val_accuracy <- trainingHistory$metrics$val_accuracy
    loss <- trainingHistory$metrics$val_loss
    out[band] <- round(val_accuracy[match(min(loss), loss)], 4)
  }

  # comparison
  selected[numberOfBands] <- names(out[match(max(out), out)])
  bands <- bands[-match(selected[numberOfBands], bands)]
}

saveRDS(selected, file.path("out", "analysis", "selectedBands.RDS"))

# when out of memory during loop
saveRDS(bands, file.path("tem", "bands.RDS"))
saveRDS(selected, file.path("tem", "selected.RDS"))
saveRDS(numberOfBands, file.path("tem", "numberOfBands.RDS"))
saveRDS(band, file.path("tem", "band.RDS"))
saveRDS(out, file.path("tem", "out.RDS"))

# quit without saving workspace image, then restart R
bands <- readRDS(file.path("tem", "bands.RDS"))
selected <- readRDS(file.path("tem", "selected.RDS"))
numberOfBands <- readRDS(file.path("tem", "numberOfBands.RDS"))
band <- readRDS(file.path("tem", "band.RDS"))
out <- readRDS(file.path("tem", "out.RDS"))

match(band, bands):length(bands)

bandsToGo <- length(selected) + 1
