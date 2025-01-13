# set parameters before plotting
mfrow_org <- par("mfrow")

newRes <- 300
f <- newRes / 72
myMar <- c(0, 0, 2, 0)

# plot
png(file.path("out", "analysis", "visualize overlap.png"), width = 480 * f, height = 250 * f, res = 72 * f, pointsize = 12)
par(mfrow = c(2, 4))

ns <- seq(0, 56, 8)
for (n in ns) {
  # read prediction raster
  overlap <- terra::rast(file.path("out", "testing", "prediction", paste0("band8 11 overlap", n), "BART_036_2019.tif"))

  # read reference raster
  r <- terra::rast(file.path("inp", "testing", "CHM", "BART_036_2019.tif"))
  overlap <- terra::crop(overlap, r)

  plot(overlap, axes = FALSE, legend = FALSE, mar = myMar, main = n, , col = .default.pal())
}
dev.off()

# recover plot parameter
par(mfrow = mfrow_org)
