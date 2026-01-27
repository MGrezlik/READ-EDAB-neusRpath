# outputting csv version of these models to use EwE keystonness calculator and compare values

# Define the data directory
data_dir <- paste0(here::here(), "/data")

# List of .rda files to export
input_files <- c(
  "GB_balanced_params.rda",
  "GOM_balanced_params.rda",
  "MAB_balanced_params.rda"
)

# Function to export model and diet components from a balanced parameter file
export_balanced_params <- function(rda_path) {
  full_path <- file.path(data_dir, rda_path)
  
  # Load into a temporary environment to avoid pollution
  temp_env <- new.env()
  load(full_path, envir = temp_env)
  
  # Get the name of the object inside the RDA file
  obj_name <- ls(envir = temp_env)
  if (length(obj_name) != 1) {
    stop(paste("Expected one object in", rda_path, "but found:", paste(obj_name, collapse = ", ")))
  }
  
  # Extract the parameters list
  params <- get(obj_name, envir = temp_env)
  
  # Verify structure
  if (!all(c("model", "diet") %in% names(params))) {
    stop(paste("File", rda_path, "does not contain expected components (model, diet)."))
  }
  
  # Convert to data frames for CSV export
  model_df <- as.data.frame(params$model)
  diet_df  <- as.data.frame(params$diet)
  
  # Derive model name (GB, GOM, MAB)
  model_name <- gsub("_balanced_params\\.rda", "", basename(rda_path))
  
  # Output CSV paths
  model_csv <- file.path(data_dir, paste0(model_name, "_model.csv"))
  diet_csv  <- file.path(data_dir, paste0(model_name, "_diet.csv"))
  
  # Write CSVs
  readr::write_csv(model_df, model_csv)
  readr::write_csv(diet_df, diet_csv)
  
  message("✅ Exported: ", model_name, " model and diet CSVs to ", data_dir)
}

# Loop through and export
for (file in input_files) {
  export_balanced_params(file)
}
