# adjust if needed
pythonDir <- "~/micromamba/envs/Chen2024/bin/python"
rootDir <- "~/Chen2024/IPO"

# set root directory
setwd(rootDir)

# create folders
foldersToCreate <- c(
  "inp", # input
  "pro", # process
  "out", # output
  "tem", # temporary
  
  "inp/training",
  "inp/training/annotations",
  "inp/training/RGB",
  "inp/training/Hyperspectral",
  "inp/training/CHM",
  
  "inp/testing",
  "inp/testing/annotations",
  "inp/testing/RGB",
  "inp/testing/Hyperspectral",
  "inp/testing/CHM",
  
  "inp/analysis",
  "inp/submission",
  
  "out/training",
  "out/training/mask",
  "out/training/annotations",
  "out/training/weight",
  "out/training/image",
  "out/training/model",
  "out/training/visualization",
  
  "out/testing",
  "out/testing/mask",
  "out/testing/image",
  "out/testing/prediction",
  "out/testing/segment",
  "out/testing/visualization",
  
  "out/analysis",
  "out/submission"
)
for (i in foldersToCreate) {
  if (!dir.exists(i)) {
    dir.create(i)
  }
}

# import functions
source("pro/function.R")

# deal with interface between python and r
assignInNamespace("is_conda_python", function(x) {
  return(FALSE)
}, ns = "reticulate")
reticulate::use_python(pythonDir)

## check if gpu is available for tensorflow
tensorflow::tf$test$is_gpu_available()

# load packages
packagesToLoad <- c(
  "terra",
  "tensorflow",
  "keras",
  "reticulate",
  "sf",
  "rsample",
  "NeonTreeEvaluation",
  "lidR",
  "parallel"
)
for (i in packagesToLoad) {
  library(i, character.only = TRUE)
}

# specify package settings
## set terra temp path
terra::terraOptions(tempdir = "tem")

## supress .aux files in terra::writeRaster
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

## use mixed precision to save memory if possible
keras::keras$mixed_precision$set_global_policy("mixed_float16")

# set common variables
## image size in pixels
targetShape <- c(128, 128)

## batch size
batch_size <- 1
