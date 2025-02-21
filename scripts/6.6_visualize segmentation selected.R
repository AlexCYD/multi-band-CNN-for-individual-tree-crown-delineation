# for paper figure
mfrow_org <- par("mfrow")

allFiles <- c("SJER_056_2018.tif", "DSNY_002_2018.tif", "TEAK_052_2018.tif", "ABBY_076_2019.tif")
labels <- c("a) SJER", "b) DSNY", "c) TEAK", "d) ABBY")

newRes <- 300
f <- newRes / 72

png(file.path(file.path("out", "analysis", "delineation.png")), width = 400 * f, height = 200 * f, res = 72 * f, pointsize = 12)
par(mfcol = c(2, 4), oma = c(0, 1, 1, 0), family = "serif")
for (i in allFiles) {
  a <- tools::file_path_sans_ext(i)

  # plot RGB
  dopPath <- file.path("inp", "testing", "RGB", i)
  dop <- terra::rast(dopPath)
  plotRGB(dop, mar = c(0.5, 0, 0, 0.5))

  # add ground truth
  ## convert xml to sf
  im <- raster::stack(dopPath)
  xmlPath <- file.path("inp", "testing", "annotations", paste0(a, ".xml"))
  dat <- dplyr::bind_rows(parallel::mclapply(xmlPath, NeonTreeEvaluation::xml_parse))
  mask <- NeonTreeEvaluation::boxes_to_spatial_polygons(dat, im)
  ## plot
  plot(mask[, "crown_id"], add = TRUE, col = NA)

  # plot prediction
  path <- file.path("out", "testing", "prediction", "band8 11 overlap24", i)
  prediction <- terra::rast(path)
  prediction <- terra::crop(prediction, dop)
  plot(prediction, axes = FALSE, mar = c(0, 0, 0.5, 0.5), legend = FALSE, col = .default.pal())

  # add segmented result
  path <- file.path("out", "testing", "segment", "band8 11 overlap24 segment76", paste0(a, ".gpkg"))
  if (file.exists(path)) {
    v <- terra::vect(path)
    plot(v, add = TRUE, border = "yellow")
  }
  
  # add title
  j <- match(i, allFiles)
  mtext(labels[j], side = 3, outer = TRUE, at = (j - 0.5) / 4, cex = 0.5)
}

# add image type
mtext("RGB + annotation", side = 2, outer = TRUE, at = 0.75, cex = 0.5)
mtext("Prediction + delineation", side = 2, outer = TRUE, at = 0.25, cex = 0.5)

dev.off()

par(mfrow = mfrow_org)
