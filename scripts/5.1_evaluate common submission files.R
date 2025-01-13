# are plots of submissions all the same
submissions <- list.files(file.path("out", "submission"), full.name = TRUE)
for (submission in submissions) {
  f <- tools::file_path_sans_ext(basename(submission))
  path <- file.path("out", "submission", paste0(f, ".RDS"))
  results <- readRDS(path)
  assign(substr(f, 13, nchar(f)), results)
}

table(Dalponte2016$plot_level$plot_name == Li2012$plot_level$plot_name)
table(Dalponte2016$plot_level$plot_name == Silva2016$plot_level$plot_name)
table(Dalponte2016$plot_level$plot_name == Weinstein_unpublished$plot_level$plot_name)
table(Dalponte2016$plot_level$plot_name == Weinstein2019$plot_level$plot_name)
# TRUE
#   63

# read my results
path <- file.path("out", "testing", "segment", "band8 11 overlap24 segment76", "imageCrowns_submissionAll.RDS")
This_study <- readRDS(path)

# are plots in submissions included in my results (yes, but just make sure)
table(Dalponte2016$plot_level$plot_name %in% This_study$plot_level$plot_name)
# TRUE
#   63

# which sites are they
a <- strsplit(Dalponte2016$plot_level$plot_name, "_", fixed = TRUE)
a <- do.call(rbind, a)
a <- a[, 2]
a <- data.frame(a)
names(a) <- c("site")
dfSites <- dplyr::count(a, site)
names(dfSites) <- c("site", "files.compare")
dfSites
#   site files
# 1 SJER    30
# 2 TEAK    33

write.csv(dfSites, file.path("out", "analysis", "dfCompare.csv"), row.names = FALSE)

# get only common plots for my submission
## overlap
n <- 24
## combination of segment parameters
j <- 76

# get common file names of the evaluation data
commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))
commonPredict <- commonPredict[commonPredict %in% Dalponte2016$plot_level$plot_name]

# generate submission file
segmentDir <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j))
inpPaths <- list.files(segmentDir, pattern = ".gpkg", full.names = TRUE)
inpPaths <- inpPaths[tools::file_path_sans_ext(basename(inpPaths)) %in% commonPredict]
generateSubmission(inpPaths, segmentDir, shouldHave = commonPredict, outName = "submissionCommon")

# evaluate submission file
segmentDir <- file.path("out", "testing", "segment", paste0("band8 11 overlap", n, " segment", j))
submission <- file.path(segmentDir, "submissionCommon.csv")
evaluateSubmission(submission, segmentDir)
