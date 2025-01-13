# overlap
n <- 24

# combination of segment parameters
j <- 76

# get common file names of the testing data
commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))

# randomly choose 19
set.seed(42)
sampling <- sample(seq_along(commonPredict), 19)

# segment
params <- expand.grid(ws = c(5, 10, 15), hmin = c(0.4, 0.5, 0.6), max_cr = c(30, 50, 70), minCrownsArea = c(0, 1.5, 3, 4.5))

for (i in commonPredict[-sampling]) {
  inpPath <- file.path("out", "testing", "prediction", paste0("band8 11 overlap", n), paste0(i, ".tif"))
  referencePath <- file.path("inp", "testing", "RGB", paste0(i, ".tif"))

  outPath <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j), paste0(i, ".gpkg"))
  segment(
    inpPath, outPath,
    referencePath = referencePath,
    ws = params$ws[j], hmin = params$hmin[j], th_tree = params$hmin[j], max_cr = params$max_cr[j], minCrownsArea = params$minCrownsArea[j]
  )
}

# generate submission file
segmentDir <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j))
inpPaths <- list.files(segmentDir, pattern = ".gpkg", full.names = TRUE)
generateSubmission(inpPaths, segmentDir, shouldHave = commonPredict, outName = "submissionAll")

# evaluate submission file
segmentDir <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j))
submission <- file.path(segmentDir, "submissionAll.csv")
evaluateSubmission(submission, segmentDir)

# view results
path <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j), "imageCrowns_submissionAll.RDS")
results <- readRDS(path)
precision <- results$overall$precision
recall <- results$overall$recall
F1score <- 2 * precision * recall / (precision + recall)

F1score
# [1] 0.3905923
precision
# [1] 0.3633763
recall
# [1] 0.4222151
