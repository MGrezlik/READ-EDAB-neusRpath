#' Calculate keystoneness indices
#'
#' Computes Keystone index #1, #2, #3 and Relative Total Impact (raw and normalized)
#' for a Northeast U.S. Rpath model, excluding detritus and fleet groups.
#'
#' @param epu Character. Ecosystem production unit: "GB", "GOM", or "MAB".
#' @return A data.frame with ecological groups, biomass, proportion biomass,
#' raw and normalized Relative Total Impact, and keystone indices.
#' @examples
#' keystoneness_all("GB")
#' @export
keystoneness_all <- function(epu = c("GB", "GOM", "MAB")) {
  
  epu <- match.arg(epu)
  
  # --- STEP 1: load model and balanced parameters ---
  data(list = epu)
  model <- get(epu)
  
  params_name <- paste0(epu, "_balanced_params")
  data(list = params_name)
  params <- get(params_name)
  
  # --- STEP 2: compute MTI matrix using local MTI function ---
  mti_mat <- MTI(Rpath = model, Rpath.params = params)
  
  # --- STEP 3: filter for living ecological groups only ---
  eco_rows <- which(params$model$Type < 2)  # living groups only
  mti_mat <- mti_mat[eco_rows, eco_rows, drop = FALSE]
  
  biomass <- model$Biomass[eco_rows]
  prop_biomass <- biomass / sum(biomass, na.rm = TRUE)
  
  # Group names
  group_names <- model$Group[eco_rows]
  if (!is.null(rownames(model$DC))) {
    group_names <- rownames(model$DC)[eco_rows]
  }
  
  # --- STEP 4: Relative Total Impact (RTI) ---
  total_effect <- rowSums(abs(mti_mat), na.rm = TRUE)
  self_effect  <- diag(abs(mti_mat))
  rti_raw      <- total_effect - self_effect
  rti_norm     <- rti_raw / max(rti_raw, na.rm = TRUE)
  
  # --- STEP 5: Keystone index #1 (Libralato 2006) ---
  KS1 <- log(rti_raw * (1 - prop_biomass))
  
  # --- STEP 6: Keystone index #2 (Valls 2015) ---
  predator_impacts <- rowSums(abs(mti_mat * (mti_mat < 0)), na.rm = TRUE)
  prey_impacts     <- rowSums(abs(mti_mat * (mti_mat > 0)), na.rm = TRUE)
  KS2 <- (predator_impacts + prey_impacts) / biomass
  
  # --- STEP 7: Keystone index #3 (Valls 2015) ---
  KS3 <- predator_impacts / biomass
  
  # --- STEP 8: tidy output ---
  results <- data.frame(
    group = group_names,
    biomass = biomass,
    prop_biomass = prop_biomass,
    RTI_raw = rti_raw,
    RTI_normalized = rti_norm,
    `Keystone index #1` = KS1,
    `Keystone index #2` = KS2,
    `Keystone index #3` = KS3,
    stringsAsFactors = FALSE
  )
  
  # Attach EPU label as attribute
  attr(results, "epu") <- switch(epu,
                                 GB  = "Georges Bank",
                                 GOM = "Gulf of Maine",
                                 MAB = "Mid Atlantic Bight"
  )
  
  return(results)
}
