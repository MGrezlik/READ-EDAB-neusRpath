library(neusRpath)

out_dir <- file.path("inst", "extdata", "keystoneness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (epu in c("GB", "GOM", "MAB")) {
  ks <- keystoneness_all(epu)
  
  write.csv(
    ks,
    file = file.path(out_dir, paste0(epu, "_keystoneness.csv")),
    row.names = FALSE
  )
}
