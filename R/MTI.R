#' Calculate Mixed Trophic Impact (MTI)
#'
#' @param Rpath Rpath object
#' @param Rpath.params Rpath parameter object
#' @param increase Logical; if TRUE, assumes marginal increase
#'
#' @return MTI matrix
#' @export
MTI <- function(Rpath, Rpath.params, increase = TRUE) {
  # --- Step 0: prepare ---
  x <- Rpath.params
  y <- Rpath
  
  all_groups <- x$model$Group
  n_groups <- length(all_groups)
  
  # --- Step 1: Build diet composition (DC) matrix ---
  DC_full <- matrix(0, nrow = n_groups, ncol = n_groups,
                    dimnames = list(all_groups, all_groups))
  
  preds <- x$diet$Group
  prey_cols <- colnames(x$diet)[-1]
  
  for(i in seq_along(preds)) {
    predator <- preds[i]
    if(!predator %in% all_groups) next
    
    prey_vals <- as.numeric(dplyr::select(x$diet[i, ], dplyr::all_of(prey_cols)))
    
    # keep only prey that exist in all_groups
    valid_idx <- which(prey_cols %in% all_groups)
    valid_preys <- prey_cols[valid_idx]
    valid_vals <- prey_vals[valid_idx]
    
    # replace NA with 0
    valid_vals[is.na(valid_vals)] <- 0
    
    DC_full[predator, valid_preys] <- valid_vals
  }
  
  # --- Step 2: Compute fraction of prey consumed by predator (FC) ---
  bio <- y$BB        # biomass
  QB  <- y$QB        # consumption/biomass
  BQB <- bio * QB    # production consumed
  
  # Tij = consumption of prey j by predator i
  Tij <- DC_full * matrix(BQB, nrow = n_groups, ncol = n_groups, byrow = TRUE)
  
  Tim <- rowSums(Tij, na.rm = TRUE)  # handle NAs
  FCij <- matrix(0, nrow = n_groups, ncol = n_groups)
  
  for(i in 1:n_groups) {
    if(!is.na(Tim[i]) && Tim[i] > 0) {
      FCij[i, ] <- Tij[i, ] / Tim[i]
    }
  }
  
  FCji <- t(FCij)
  
  # --- Step 3: Net impact ---
  net_impact <- DC_full - FCji
  
  # --- Step 4: MTI ---
  MTI <- MASS::ginv(diag(n_groups) - net_impact) - diag(n_groups)
  
  return(MTI)
}
