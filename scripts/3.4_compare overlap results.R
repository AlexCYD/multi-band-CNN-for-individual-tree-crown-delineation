# collect accuracy from overlap models
out <- c()
ns <- seq(0, 56, 8)
for (n in ns) {
  path <- file.path("out", "testing", "prediction", paste0("band8 11 overlap", n), "BART_036_2019.RDS")
  out <- c(out, round(readRDS(path)["accuracy"], 4))
}
out
# accuracy accuracy accuracy accuracy accuracy accuracy accuracy accuracy 
#   0.7208   0.7367   0.7412   0.7525   0.7553   0.7508   0.7598   0.7645

# plot
newRes <- 300
f <- newRes / 72

png(file.path("out", "analysis", "compare overlap.png"), width = 400 * f, height = 200 * f, res = 72 * f, pointsize = 12)
par(mar = c(4, 4, 0, 0) + 0.1, family = "serif")
plot(ns, out, xlab = "Overlap distance (pixel)", ylab = "Test accuracy", type = "b", xaxt = "n")
axis(1, at = ns)
dev.off()
