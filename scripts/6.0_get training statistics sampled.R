# How many tiles are there per site?
# only sample some due to memory issue later in training
maskDir <- file.path("out", "training", "mask")
set.seed(42)
sampling <- sample(seq_along(list.files(maskDir)), 100)
a <- list.files(maskDir)
a <- a[sampling]
a <- strsplit(a, "_", fixed = TRUE)
a <- do.call(rbind, a)
a <- a[, c(1:3)]
a <- data.frame(a)
names(a) <- c("year", "site", "plot")
dfTiles <- dplyr::count(a, site)
names(dfTiles) <- c("site", "tiles.Sample")
dfTiles

write.csv(dfTiles, file.path("out", "analysis", "dfTrainSample.csv"), row.names = FALSE)
