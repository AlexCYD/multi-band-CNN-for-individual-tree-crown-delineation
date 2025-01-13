# How many files are there per site?
# sample those for tuning
commonPredict <- tools::file_path_sans_ext(list.files(file.path("inp", "testing", "RGB")))
set.seed(42)
sampling <- sample(seq_along(commonPredict), 19)
a <- sort(commonPredict[sampling])

a1 <- a[1:4]
a2 <- a[5:19]

a1 <- strsplit(a1, "_", fixed = TRUE)
a1 <- do.call(rbind, a1)
a1 <- a1[, c(1:3)]
a1 <- data.frame(a1)
names(a1) <- c("year", "site", "plot")

a2 <- strsplit(a2, "[._]")
a2 <- do.call(rbind, a2)
a2 <- a2[, c(1:3)]
a2 <- data.frame(a2)
names(a2) <- c("site", "plot", "year")

a <- rbind(a1, a2)

dfFiles <- dplyr::count(a, site)
names(dfFiles) <- c("site", "files.sample")
dfFiles

write.csv(dfFiles, file.path("out", "analysis", "dfEvalSample.csv"), row.names = FALSE)
