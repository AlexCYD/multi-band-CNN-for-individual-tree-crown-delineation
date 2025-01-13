# overlap
n <- 24
overlapStride <- c(n, n)

# get common file names
commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))

# set parameters for myPredict
modelDir <- file.path("out", "training", "model", "band8 11")
outDir <- file.path("out", "testing", "prediction", paste0("band8 11 overlap", n))
model <- unet(c(targetShape, 8), weighted = FALSE, num_layers = 4)
loss <- "binary_crossentropy"
metrics <- "accuracy"

# get file names without tile numbering
imageDir <- file.path("out", "testing", "image", paste0("overlap", n))
maskDir <- file.path("out", "testing", "mask", paste0("overlap", n))
files <- data.frame(
  imagePaths = list.files(imageDir, full.names = TRUE),
  maskPaths = list.files(maskDir, full.names = TRUE)
)
filesWithoutNumbering <- files
filesWithoutNumbering$imagePaths <- gsub("_[0-9]+(?=\\.tif$)", "", basename(files$imagePaths), perl = TRUE)

# generate test and prediction per testing file
for (i in commonPredict[-66]) {
  selected <- files[filesWithoutNumbering$imagePaths == paste0(i, ".tif"), ]
  imagePaths <- selected$imagePaths
  maskPaths <- selected$maskPaths
  optimizer <- keras::optimizer_sgd(learning_rate = 0.001, momentum = 0.99)
  myPredict(
    imagePaths, maskPaths, modelDir,
    model, optimizer, loss, metrics,
    batch_size, outDir, overlapStride,
    doEvaluation = FALSE
  )
}
