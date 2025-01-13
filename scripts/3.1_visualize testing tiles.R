# overlap pixel
ns <- seq(0, 56, 8)

# set parameters before plotting
mfrow_org <- par("mfrow")

newRes <- 300
f <- newRes / 72
myMar <- c(0, 0, 2, 0)

for (n in ns) {
  imageDir <- file.path("out", "testing", "image", paste0("overlap", n))
  allFiles <- list.files(imageDir)

  # get common file names of the testing data
  commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))
  allFiles <- allFiles[grepl(commonPredict[66], allFiles, fixed = TRUE)]

  for (i in allFiles) {
    outDir <- file.path("out", "testing", "visualization", paste0("image_overlap", n))
    if (!dir.exists(outDir)) {
      dir.create(outDir)
    }

    # read file
    imagePath <- file.path(imageDir, i)
    imager <- terra::rast(imagePath)

    # plot every band
    png(file.path(outDir, paste0(tools::file_path_sans_ext(i), ".png")), width = 480 * f, height = 250 * f, res = 72 * f, pointsize = 12)
    par(mfrow = c(2, 4))
    for (j in seq_len(nlyr(imager))) {
      plot(imager[[j]], axes = FALSE, legend = FALSE, mar = myMar, main = names(imager[[j]]), col = .default.pal())
    }
    dev.off()
  }
}

# recover plot parameter
par(mfrow = mfrow_org)
