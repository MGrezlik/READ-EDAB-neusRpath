# outputting csv version of these models to use EwE keystonness calculator and compare values

# Define the data directory where your .rda files live
data_dir <- paste0(here::here(),"/data")

# List of .rda files to export
input_files <- c(
  "GB_balanced_params.rda",
  "GOM_balanced_params.rda",
  "MAB_balanced_params.rda"
)

# Function to export model and diet components from a balanced parameter file
export_balanced_params <- function(rda_path) {
  full_path <- file.path(data_dir, rda_path)
  
  # Load the R object
  load(full_path)
  
  # Identify the object name (e.g., "GB_balanced_params")
  obj_name <- ls()[1]
  params <- get(obj_name)
  
  # Ensure expected structure exists
  if (!all(c("model", "diet") %in% names(params))) {
    stop(paste("File", rda_path, "does not contain expected components."))
  }
  
  # Convert to data frames for CSV export
  model_df <- as.data.frame(params$model)
  diet_df  <- as.data.frame(params$diet)
  
  # Derive model name (GB, GOM, MAB)
  model_name <- gsub("_balanced_params\\.rda", "", basename(rda_path))
  
  # Output CSV paths
  model_csv <- file.path(data_dir, paste0(model_name, "_model.csv"))
  diet_csv  <- file.path(data_dir, paste0(model_name, "_diet.csv"))
  
  # Write CSV files
  readr::write_csv(model_df, model_csv)
  readr::write_csv(diet_df, diet_csv)
  
  message("✅ Exported: ", model_name, " model and diet CSVs to ", data_dir)
}

# Loop through each file and export
for (file in input_files) {
  export_balanced_params(file)
}
