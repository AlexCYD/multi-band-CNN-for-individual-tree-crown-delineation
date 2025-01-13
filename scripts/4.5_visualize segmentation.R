outDir <- file.path("out", "testing", "visualization", "band8 11 overlap24 segment76")
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

mfrow_org <- par("mfrow")

allFiles <- list.files(file.path("inp", "testing", "RGB"))

for (i in allFiles) {
  a <- tools::file_path_sans_ext(i)

  newRes <- 300
  f <- newRes / 72

  png(file.path(outDir, paste0(a, ".png")), width = 245 * f, height = 400 * f, res = 72 * f, pointsize = 12)
  par(mfrow = c(2, 1))

  # plot RGB
  dopPath <- file.path("inp", "testing", "RGB", i)
  dop <- terra::rast(dopPath)
  plotRGB(dop, mar = c(0.5, 0, 0, 3))

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
  plot(prediction, axes = FALSE, mar = c(0, 0, 0.5, 3))

  # add segmented result
  path <- file.path("out", "testing", "segment", "band8 11 overlap24 segment76", paste0(a, ".gpkg"))
  if (file.exists(path)) {
    v <- terra::vect(path)
    plot(v, add = TRUE, border = "yellow")
  }
  dev.off()
}

par(mfrow = mfrow_org)
