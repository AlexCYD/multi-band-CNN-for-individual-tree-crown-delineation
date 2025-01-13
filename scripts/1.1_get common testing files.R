# adjust if needed
testingDir <- file.path("..", "NeonTreetesting")

# define paths
annotationsDir <- file.path(testingDir, "annotations")
rgbDir <- file.path(testingDir, "testing", "RGB")
hyperspectralDir <- file.path(testingDir, "testing", "Hyperspectral")
chmDir <- file.path(testingDir, "testing", "CHM")

common <- getCommonFiles(annotationsDir, rgbDir, hyperspectralDir, chmDir)

# annotations
files_to_copy <- file.path(annotationsDir, paste0(common, ".xml"))
outDir <- file.path("inp", "testing", "annotations")

file.copy(files_to_copy, outDir, overwrite = TRUE)

# RGB
files_to_copy <- file.path(rgbDir, paste0(common, ".tif"))
outDir <- file.path("inp", "testing", "RGB")

file.copy(files_to_copy, outDir, overwrite = TRUE)

# hyperspectral
files_to_copy <- file.path(hyperspectralDir, paste0(common, "_hyperspectral.tif"))
outDir <- file.path("inp", "testing", "Hyperspectral", paste0(common, ".tif"))

file.copy(files_to_copy, outDir, overwrite = TRUE)

# CHM
files_to_copy <- file.path(chmDir, paste0(common, "_CHM.tif"))
outDir <- file.path("inp", "testing", "CHM", paste0(common, ".tif"))

file.copy(files_to_copy, outDir, overwrite = TRUE)
