#' Calculate Fishable Production with Automated Area Calculation
#'
#' @param model A balanced Rpath model object.
#' @param ecosystem_name Character. The name of the ecosystem. 
#'        Accepted values: "Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight".
#' @param base_year Numeric. Base year of the model.
#' @param citation Character. Citation for the model.
#' @param exclude_groups Character vector. Names of groups to exclude.
#' @param area_km2 Numeric. Optional. If provided, overrides the shapefile calculation.
#' @param shapefile_path Character. Path to the EPU shapefile.
#' @param fleet_name Character. Optional name of the dominant fleet.
#' @param kappa Numeric. Harvesting factor. Default is 0.175.
#' @export
calc_fishable_production <- function(model, ecosystem_name, base_year, 
                                     citation, exclude_groups = c(), 
                                     area_km2 = NULL,
                                     shapefile_path = "data-raw/EPU_extended.shp",
                                     fleet_name = NULL, kappa = 0.175) {
  
  # 1. Automated Area Calculation (if area_km2 is not provided)
  if (is.null(area_km2)) {
    if (!file.exists(shapefile_path)) {
      stop("Shapefile not found at ", shapefile_path, ". Please provide area_km2 manually or check the path.")
    }
    
    # Map user-friendly names to the EPU codes in your shapefile
    epu_map <- c("Georges Bank" = "GB", 
                 "Gulf of Maine" = "GOM", 
                 "Mid-Atlantic Bight" = "MAB")
    
    target_epu <- epu_map[ecosystem_name]
    
    if (is.na(target_epu)) {
      stop("ecosystem_name must be 'Georges Bank', 'Gulf of Maine', or 'Mid-Atlantic Bight' to use automated area.")
    }
    
    # Temporarily disable strict s2 spherical geometry to bypass vertex errors
    current_s2 <- sf::sf_use_s2()
    suppressMessages(sf::sf_use_s2(FALSE))
    
    # Read, fix geometry, and calculate area
    epu_sf <- sf::st_read(shapefile_path, quiet = TRUE)
    
    area_km2 <- epu_sf %>%
      dplyr::filter(EPU == target_epu) %>%
      sf::st_make_valid() %>% # Repairs self-intersections
      sf::st_area() %>%
      as.numeric() / 1000000 # Convert m2 to km2
    
    # Restore s2 to its original state so we don't mess up the user's environment
    suppressMessages(sf::sf_use_s2(current_s2))
    
    message(paste("Automated Area Lookup:", round(area_km2, 2), "km2 for", ecosystem_name))
  }
  
  # 2. Extract parameters from Rpath model
  groups  <- model$Group
  TL      <- model$TL
  Biomass <- model$Biomass
  PB      <- model$PB
  catch_mat <- model$Landings + model$Discards
  
  # 3. Determine dominant fleet
  fleet_totals <- colSums(catch_mat, na.rm = TRUE)
  if (is.null(fleet_name)) {
    fleet_name <- names(fleet_totals)[which.max(fleet_totals)]
  }
  
  # 4. Calculate Fleet Trophic Level
  fleet_catch <- catch_mat[, fleet_name]
  fleet_TL <- sum(fleet_catch * TL, na.rm = TRUE) / sum(fleet_catch, na.rm = TRUE)
  
  # 5. Set TL Threshold & 90% Catch Coverage Loop
  initial_tl_threshold <- fleet_TL - 1.0
  current_tl_threshold <- initial_tl_threshold
  total_catch <- sum(catch_mat, na.rm = TRUE)
  perc_catch <- 0
  
  while(perc_catch < 0.90 && current_tl_threshold > 0) {
    valid_taxa <- !(groups %in% exclude_groups) & (TL >= current_tl_threshold)
    catch_cov  <- sum(catch_mat[valid_taxa, ], na.rm = TRUE)
    perc_catch <- catch_cov / total_catch
    if(perc_catch < 0.90) current_tl_threshold <- current_tl_threshold - 0.01
  }
  
  # 6. Calculate Fishable Production
  valid_taxa_final <- !(groups %in% exclude_groups) & (TL >= current_tl_threshold)
  total_production <- sum(Biomass[valid_taxa_final] * PB[valid_taxa_final], na.rm = TRUE)
  
  fishable_prod_km2 <- total_production * kappa
  fishable_prod_total <- fishable_prod_km2 * area_km2
  
  # 7. Format Result Paragraph
  extension_text <- ""
  if (current_tl_threshold < initial_tl_threshold) {
    extension_text <- sprintf(
      " This trophic level cutoff was extended down to TL %.2f, to ensure that at least 90%% of harvest in the base model was included.", 
      current_tl_threshold
    )
  }
  
  template_text <- sprintf(
    "We applied the Fishable Production method to an Rpath food web model of %s, which represented a base year of %d. This model (%s) is intended to represent a spatial domain of %s km2. In this Rpath model, the single fleet or dominant fleet %s had a trophic level of %.2f, and hence we include production estimates of all species >= trophic level %.2f.%s We estimate Fishable production to be %.2f tons km-2, equivalent to %.2f million tons over the model spatial domain.",
    ecosystem_name, base_year, citation, format(round(area_km2), big.mark = ","),
    fleet_name, fleet_TL, initial_tl_threshold, extension_text,
    fishable_prod_km2, fishable_prod_total / 1000000
  )
  
  return(list(metrics = data.frame(Ecosystem = ecosystem_name, Area = area_km2, Fishable_Prod_km2 = fishable_prod_km2),
              report_text = template_text))
}