#' Ranked biomass diagnostics for balanced Rpath models
#'
#' Extracts living-group biomasses from a balanced Rpath model
#' and returns a ranked table (largest to smallest).
#'
#' @param rpath.obj A balanced Rpath model object
#'
#' @return A data.frame of living groups ranked by biomass
#' @export

diagnostics_balanced_biomass <- function(rpath.obj) {
  
  if (!requireNamespace("Rpath", quietly = TRUE)) {
    stop("Package 'Rpath' required.")
  }
  
  model <- as.data.frame(Rpath::write.Rpath(rpath.obj))
  
  # Living groups only (type < 2)
  living <- model[model$type < 2 & model$Biomass > 0,
                  c("Group", "TL", "Biomass")]
  
  # Order by descending biomass
  living <- living[order(living$Biomass, decreasing = TRUE), ]
  
  # Add rank column
  living$Rank <- seq_len(nrow(living))
  
  rownames(living) <- NULL
  
  return(living)
}