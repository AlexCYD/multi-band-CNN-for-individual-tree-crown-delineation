selectedBands <- readRDS(file.path("out", "analysis", "selectedBands.RDS"))

# overlap
n <- 24
overlapStride <- c(n, n)

# get common file names of the testing data
commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))

# generate testing mask
outDir <- file.path("out", "testing", "mask", paste0("overlap", n))
for (i in commonPredict[-66]) {
  xmlPath <- file.path("inp", "testing", "annotations", paste0(i, ".xml"))
  dopPath <- file.path("inp", "testing", "RGB", paste0(i, ".tif"))
  generateMask(
    xmlPath, dopPath, targetShape, outDir,
    mixedOnly = FALSE, overlapStride = overlapStride
  )
}

# generate testing image
maskDir <- file.path("out", "testing", "mask", paste0("overlap", n))
outDir <- file.path("out", "testing", "image", paste0("overlap", n))
files <- list.files(maskDir, full.names = TRUE)
filesWithoutNumbering <- gsub("_[0-9]+(?=\\.tif$)", "", basename(files), perl = TRUE)
for (i in commonPredict[-66]) {
  selected <- files[filesWithoutNumbering == paste0(i, ".tif")]

  dopPath <- file.path("inp", "testing", "RGB", paste0(i, ".tif"))
  hyperspectralPath <- file.path("inp", "testing", "Hyperspectral", paste0(i, ".tif"))
  chmPath <- file.path("inp", "testing", "CHM", paste0(i, ".tif"))

  for (maskPath in selected) {
    generateImage(
      maskPath, dopPath, hyperspectralPath, chmPath, outDir,
      overlapStride = overlapStride, useBands = selectedBands[1:8]
    )
  }
}
