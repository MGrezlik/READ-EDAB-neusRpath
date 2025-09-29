#' Calculate keystoneness for a Northeast U.S. Rpath model
#'
#' Calculates the keystoneness index (KSi) for a given ecosystem model
#' following Libralato et al. 2006:
#' KSi = log[epsilon_i * (1 - p_i)],
#' where epsilon_i is the total impact of group i (from the MTI matrix excluding self-effects)
#' and p_i is the proportion of total system biomass.
#'
#' @param epu Character. Ecosystem production unit: "GB", "GOM", or "MAB".
#' @return A data.frame with ecological groups, biomass, proportion biomass,
#' total impact, and keystoneness index.
#' @examples
#' keystoneness("GB")
#' @export
keystoneness <- function(epu = c("GB", "GOM", "MAB")) {
  
  epu <- match.arg(epu)
  
  # --- STEP 1: load model and balanced parameters ---
  data(list = paste0(epu))  # loads GB, GOM, or MAB model into environment
  model <- get(epu)
  
  # load balanced parameters
  params_name <- paste0(epu, "_balanced_params")
  data(list = params_name)
  params <- get(params_name)
  
  # --- STEP 2: Mixed trophic impact matrix ---
  mti_mat <- Rpath::MTI(Rpath = model, Rpath.params = params)
  
  # --- STEP 3: assign group names (ecological groups only) ---
  if (!is.null(rownames(model$DC))) {
    group_names <- rownames(model$DC)
  } else if (!is.null(model$Group)) {
    group_names <- model$Group
  } else {
    group_names <- paste0("Group", seq_len(nrow(mti_mat)))
  }
  
  # exclude fishing gear groups
  n_eco <- nrow(mti_mat) - model$NUM_GEARS
  mti_mat <- mti_mat[1:n_eco, 1:n_eco]
  group_names <- group_names[1:n_eco]
  
  # --- STEP 4: compute proportional biomass ---
  biomass <- model$Biomass[1:n_eco]
  prop_biomass <- biomass / sum(biomass)
  
  # --- STEP 5: compute total impacts excluding self-effects ---
  total_effect <- apply(abs(mti_mat), 1, sum)
  self_effect <- diag(abs(mti_mat))
  epsilon <- total_effect - self_effect
  
  # --- STEP 6: keystoneness (Libralato et al. 2006) ---
  KS_index <- log(epsilon * (1 - prop_biomass))
  
  # --- STEP 7: tidy output ---
  results <- data.frame(
    group = group_names,
    biomass = biomass,
    prop_biomass = prop_biomass,
    impact = epsilon,
    keystoneness = KS_index,
    stringsAsFactors = FALSE
  )
  
  return(results)
}
