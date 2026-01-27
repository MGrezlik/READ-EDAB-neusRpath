#' Calculate Ecopath-consistent keystoneness indices
#'
#' Computes Keystone indices #1, #2, #3 and Relative Total Impact
#' exactly as implemented in Ecopath (Valls et al. + legacy code).
#'
#' @param epu Character. "GB", "GOM", or "MAB"
#' @return data.frame of keystoneness metrics for living groups
#' @export
keystoneness_all <- function(epu = c("GB", "GOM", "MAB")) {
  
  epu <- match.arg(epu)
  
  ## ------------------------------------------------------------------
  ## Load model and balanced parameters
  ## ------------------------------------------------------------------
  
  data(list = epu, envir = environment())
  model <- get(epu, envir = environment())
  
  params_name <- paste0(epu, "_balanced_params")
  data(list = params_name, envir = environment())
  params <- get(params_name, envir = environment())
  
  ## ------------------------------------------------------------------
  ## Compute MTI using the updated Ecopath-consistent function
  ## ------------------------------------------------------------------
  
  mti_full <- MTI_new(Rpath = model, Rpath.params = params)
  
  ## ------------------------------------------------------------------
  ## Restrict to living groups only
  ## ------------------------------------------------------------------
  
  NUM_LIVING <- model$NUM_LIVING
  mti_living <- mti_full[1:NUM_LIVING, 1:NUM_LIVING, drop = FALSE]
  diag(mti_living) <- 0
  
  biomass <- model$Biomass[1:NUM_LIVING]
  group_names <- model$Group[1:NUM_LIVING]
  
  ## ------------------------------------------------------------------
  ## Overall effect (epsilon_i) for keystoneness
  ## ------------------------------------------------------------------
  
  epsilon_i <- sqrt(rowSums(mti_living^2))
  
  ## ------------------------------------------------------------------
  ## Relative Total Impact (Ecopath-style)
  ## ------------------------------------------------------------------
  
  # Absolute MTI values
  mti_abs <- abs(mti_living)
  
  # Total effect per group
  total_effect <- rowSums(mti_abs, na.rm = TRUE)
  
  # Self-effect (diagonal)
  self_effect <- diag(mti_abs)
  
  # Raw RTI (exclude self-effect)
  rti_raw <- total_effect - self_effect
  
  # Normalize to 0-1
  rti <- if (all(is.na(rti_raw)) || max(rti_raw, na.rm = TRUE) == 0) {
    rep(0, length(rti_raw))
  } else {
    rti_raw / max(rti_raw, na.rm = TRUE)
  }
  
  ## ------------------------------------------------------------------
  ## Keystone indices (Ecopath / Valls)
  ## ------------------------------------------------------------------
  
  eps <- .Machine$double.eps
  p_i <- biomass / sum(biomass)
  
  KS_1 <- log10(pmax(epsilon_i * (1 - p_i), eps))
  KS_2 <- log10(pmax(epsilon_i / p_i, eps))
  
  # Biomass rank (descending)
  drank <- NUM_LIVING - rank(biomass, ties.method = "average") + 1
  KS_3 <- log10(pmax(epsilon_i * drank, eps))
  
  ## ------------------------------------------------------------------
  ## Assemble output
  ## ------------------------------------------------------------------
  
  results <- data.frame(
    group                 = group_names,
    biomass               = biomass,
    prop_biomass          = p_i,
    Relative_total_impact = rti,
    Keystone_index_1      = KS_1,
    Keystone_index_2      = KS_2,
    Keystone_index_3      = KS_3,
    stringsAsFactors      = FALSE
  )
  
  attr(results, "epu") <- switch(
    epu,
    GB  = "Georges Bank",
    GOM = "Gulf of Maine",
    MAB = "Mid Atlantic Bight"
  )
  
  return(results)
}


