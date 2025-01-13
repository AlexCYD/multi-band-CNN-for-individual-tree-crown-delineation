# read training history
i <- file.path("out", "training", "model", "band8 11")
path <- file.path(i, "trainingHistory.RDS")
trainingHistory <- readRDS(path)

# plot training history
trainingHistory$params$epochs <- 25

newRes <- 300
f <- newRes / 72

png(file.path("out", "analysis", "trainingHistory final.png"), width = 400 * f, height = 400 * f, res = 72 * f, pointsize = 12, family = "serif")
    plot(trainingHistory, method = "ggplot2", theme_bw = TRUE)
dev.off()
