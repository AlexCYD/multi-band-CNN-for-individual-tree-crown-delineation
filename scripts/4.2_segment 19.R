# overlap
n <- 24

# get common file names of the evaluation data
commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))

# randomly choose 19
set.seed(42)
sampling <- sample(seq_along(commonPredict), 19)
commonPredict <- commonPredict[sampling]

# segment with grid search
params <- expand.grid(ws = c(5, 10, 15), hmin = c(0.4, 0.5, 0.6), max_cr = c(30, 50, 70), minCrownsArea = c(0, 1.5, 3, 4.5))

for (i in commonPredict) {
  inpPath <- file.path("out", "testing", "prediction", paste0("band8 11 overlap", n), paste0(i, ".tif"))
  referencePath <- file.path("inp", "testing", "RGB", paste0(i, ".tif"))

  for (j in seq_len(nrow(params))) {
    outPath <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j), paste0(i, ".gpkg"))
    segment(
      inpPath, outPath,
      referencePath = referencePath,
      ws = params$ws[j], hmin = params$hmin[j], th_tree = params$hmin[j], max_cr = params$max_cr[j], minCrownsArea = params$minCrownsArea[j]
    )
  }
}

# generate submission file
for (j in seq_len(nrow(params))) {
  segmentDir <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j))
  inpPaths <- list.files(segmentDir, pattern = ".gpkg", full.names = TRUE)
  inpPaths <- inpPaths[tools::file_path_sans_ext(basename(inpPaths)) %in% commonPredict]
  generateSubmission(inpPaths, segmentDir, shouldHave = commonPredict, outName = "submission")
}

# evaluate submission file
for (j in seq_len(nrow(params))) {
  segmentDir <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j))
  submission <- file.path(segmentDir, "submission.csv")
  evaluateSubmission(submission, segmentDir)
}
