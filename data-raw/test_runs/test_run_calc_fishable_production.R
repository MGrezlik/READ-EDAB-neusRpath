# GB ----------

# 1. Load required tools
library(sf) # Required for the automated area calculation
library(dplyr)

# 2. Load the function and the data
# If you are actively building the package, devtools::load_all() is best.
# Alternatively, you can just source the function directly:
source("R/calc_fishable_production.R")

# Load the Georges Bank balanced model object
load("data/GB.rda")

# 3. Define the groups to exclude
exclude_list <- c(
  "SeaBirds", "Pinnipeds", "BaleenWhales", "Odontocetes", 
  "Bacteria", "Detritus", "Phytoplankton", 
  "LgCopepods", "SmCopepods", "Microzooplankton", 
  "GelZooplankton", "Krill", "Micronekton", 
  "Macrobenthos", "Megabenthos", "Fauna"
)

# 4. Run the function
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

gb_test_run <- calc_fishable_production(
  model = GB, # The object loaded from GB.rda
  ecosystem_name = "Georges Bank",
  base_year = 1985,
  citation = "Weisberg et al. in review",
  exclude_groups = exclude_list
)

Sys.unsetenv("SHAPE_RESTORE_SHX")

# 5. View the outputs
print(gb_test_run$metrics)

cat("\n--- Generated Template Paragraph ---\n")
cat(gb_test_run$report_text, "\n")

# GOM --------

# Load the Gulf of Maine balanced model object
load("data/GOM.rda")

# 4. Run the function
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

gom_test_run <- calc_fishable_production(
  model = GOM,
  ecosystem_name = "Gulf of Maine",
  base_year = 1985,
  citation = "Weisberg et al. in review",
  exclude_groups = exclude_list
)

Sys.unsetenv("SHAPE_RESTORE_SHX")

# 5. View the outputs
print(gom_test_run$metrics)

cat("\n--- Generated Template Paragraph ---\n")
cat(gom_test_run$report_text, "\n")

# MAB ---------

# Load the Gulf of Maine balanced model object
load("data/MAB.rda")


# 4. Run the function
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

mab_test_run <- calc_fishable_production(
  model = MAB,
  ecosystem_name = "Mid-Atlantic Bight",
  base_year = 1985,
  citation = "Weisberg et al. in review",
  exclude_groups = exclude_list
)

Sys.unsetenv("SHAPE_RESTORE_SHX")

# 5. View the outputs
print(mab_test_run$metrics)

cat("\n--- Generated Template Paragraph ---\n")
cat(mab_test_run$report_text, "\n")
