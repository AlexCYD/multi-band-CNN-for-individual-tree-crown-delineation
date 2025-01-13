# get all filenames
maskDir <- file.path("out", "training", "mask")
set.seed(42)
sampling <- sample(seq_along(list.files(maskDir)), 100)
allFiles <- list.files(maskDir)
allFiles <- allFiles[sampling]

# set parameters before plotting
mfrow_org <- par("mfrow")
newRes <- 300
f <- newRes / 72
myMar <- c(0, 0, 2, 0)

for (i in allFiles) {
  # read files
  maskPath <- file.path("out", "training", "mask", i)
  imagePath <- file.path("out", "training", "image", i)
  weightPath <- file.path("out", "training", "weight", i)
  imager <- terra::rast(imagePath)
  maskr <- terra::rast(maskPath)
  weightr <- terra::rast(weightPath)

  # plot
  filename <- file.path("out", "training", "visualization", paste0(tools::file_path_sans_ext(i), ".png"))
  png(filename, width = 480 * f, height = 600 * f, res = 72 * f, pointsize = 12)
  par(mfrow = c(5, 4))
  for (j in seq_len(nlyr(imager))) {
    plot(imager[[j]], axes = FALSE, legend = FALSE, mar = myMar, main = names(imager[[j]]), col = .default.pal())
  }
  plot.new()
  plotRGB(imager[[1:3]] * 255, mar = myMar, main = "RGB")
  plotRGB(imager[[4:6]] * 255, mar = myMar, main = "11/55/113")
  plot(maskr, axes = FALSE, legend = FALSE, mar = myMar, main = "target", col = .default.pal())
  plot(weightr, axes = FALSE, legend = FALSE, mar = myMar, main = "weight", col = .default.pal())
  dev.off()
}

# recover plot parameter
par(mfrow = mfrow_org)
