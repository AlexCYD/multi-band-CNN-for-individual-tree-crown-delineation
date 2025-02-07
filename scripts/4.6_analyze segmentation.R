# show evaluation result
path <- file.path("out", "testing", "segment", "band8 11 overlap24 segment76", "imageCrowns_submissionAll.RDS")
results <- readRDS(path)
results
# $overall
# A tibble: 1 × 2
#   precision recall
#       <dbl>  <dbl>
# 1     0.363  0.422

# $by_site
# A tibble: 20 × 3
# Groups:   Site [20]
#    Site  recall precision
#    <chr>  <dbl>     <dbl>
#  1 ABBY   0.013     0.028
#  2 BART   0.237     0.31 
#  3 BLAN   0.274     0.299
#  4 BONA   0.09      0.161
#  5 CLBJ   0.155     0.281
#  6 DELA   0.23      0.4  
#  7 DSNY   0.322     0.475
#  8 HARV   0.146     0.212
#  9 JERC   0.376     0.392
# 10 LENO   0.293     0.306
# 11 MLBS   0.21      0.341
# 12 NIWO   0.01      0.035
# 13 OSBS   0.207     0.333
# 14 SCBI   0.233     0.254
# 15 SERC   0.202     0.264
# 16 SJER   0.632     0.388
# 17 SOAP   0.105     0.25 
# 18 TALL   0.419     0.433
# 19 TEAK   0.221     0.333
# 20 WREF   0.281     0.318

precision <- results$overall$precision
precision
# [1] 0.3633763
recall <- results$overall$recall
recall
# [1] 0.4222151
F1score <- 2 * precision * recall / (precision + recall)
F1score
# [1] 0.3905923

# plot precision vs recall 
library(ggplot2)
library(ggrepel)

map_plot <- ggplot() +
  geom_point(data = results$by_site, aes(x = precision, y = recall), shape = 1, size = 1.5) + # Points
  ylim(0, 0.65) +
  xlim(0, 0.65) +
  geom_text_repel(data = results$by_site, vjust = 0.2, aes(x = precision, y = recall, label = Site), size = 3.8, family = "Times New Roman") + # Labels
  labs(title = NULL, x = "Precision", y = "Recall") + 
  theme_minimal() + # A clean theme for the map
  theme(
    panel.grid = element_blank(), # remove grid
    panel.border = element_rect(color = "black", fill = NA), # Add black frame
    axis.ticks = element_line(), # add tick
    axis.ticks.length = unit(0.15, "cm"), # outer tick
    axis.text = element_text(size = 11, color = "black"), # axis text
    text = element_text(family = "Times New Roman", color = "black"),
    plot.margin = unit(c(0, 0, 0, 0), "cm") # remove margin
  )

f <- file.path("out", "analysis", "evaluation per site.png")
ggsave(f, plot = map_plot, dpi = 300, width = 5, height = 5)
