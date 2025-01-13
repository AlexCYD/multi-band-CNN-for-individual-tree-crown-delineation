# calculate F1-score
path <- file.path("out", "testing", "segment", "band8 11 overlap24 segment76", "imageCrowns_submissionAll.RDS")
results <- readRDS(path)

results$by_site$F1score <- 2 * results$by_site$precision * results$by_site$recall / (results$by_site$precision + results$by_site$recall)

# calculate site density
dfEval <- read.csv(file.path("out", "analysis", "dfEval.csv"))
dfEval$density <- dfEval$crowns / dfEval$files / (0.4 * 0.4) # crown per km^2

# merge data
df <- merge(dfEval, results$by_site, by.x = "site", by.y = "Site")

# plot F1 vs tree crown density
newRes <- 300
f <- newRes / 72

png(file.path("out", "analysis", "F1 vs crown density.png"), width = 400 * f, height = 400 * f, res = 72 * f, pointsize = 12)
par(mar = c(4, 4, 0, 0) + 0.1, family = "serif")
plot(df$density, df$F1score, type = "n",
  xlab = expression(Cronws~per~km^2), ylab = "F1-score",
  ylim = c(0, 0.5), xlim = c(0, 1000)
)
text(df$density, df$F1score, df$site)
dev.off()
