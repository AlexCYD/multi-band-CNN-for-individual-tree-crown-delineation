# get common file names of the training data
commonTrain <- tools::file_path_sans_ext(list.files(file.path("inp", "training", "RGB")))

# generate training mask
outDir <- file.path("out", "training", "mask")
gpkgDir <- file.path("out", "training", "annotations")
for (i in commonTrain) {
  xmlPath <- file.path("inp", "training", "annotations", paste0(i, ".xml"))
  dopPath <- file.path("inp", "training", "RGB", paste0(i, ".tif"))
  generateMask(
    xmlPath, dopPath, targetShape, outDir,
    gpkgDir = gpkgDir, mixedOnly = TRUE
  )
}

length(list.files(outDir)) # 9061

# only sample some due to memory issue later in training
outDir <- file.path("out", "training", "mask")
set.seed(42)
sampling <- sample(seq_along(list.files(outDir)), 100)

# generate training weight
maskDir <- file.path("out", "training", "mask")
outDir <- file.path("out", "training", "weight")
files <- list.files(maskDir, full.names = TRUE)
for (i in sampling) {
  maskPath <- files[i]

  a <- tools::file_path_sans_ext(gsub("_[0-9]+(?=\\.tif$)", "", basename(maskPath), perl = TRUE))
  gpkgPath <- file.path("out", "training", "annotations", paste0(a, ".gpkg"))

  generateWeight(maskPath, gpkgPath, outDir)
}

# generate training image
maskDir <- file.path("out", "training", "mask")
outDir <- file.path("out", "training", "image")
files <- list.files(maskDir, full.names = TRUE)
for (i in sampling) {
  maskPath <- files[i]

  a <- tools::file_path_sans_ext(gsub("_[0-9]+(?=\\.tif$)", "", basename(maskPath), perl = TRUE))
  dopPath <- file.path("inp", "training", "RGB", paste0(a, ".tif"))
  hyperspectralPath <- file.path("inp", "training", "Hyperspectral", paste0(a, ".tif"))
  chmPath <- file.path("inp", "training", "CHM", paste0(a, ".tif"))

  generateImage(maskPath, dopPath, hyperspectralPath, chmPath, outDir)
}
