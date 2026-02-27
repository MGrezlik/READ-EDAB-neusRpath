#' PREBAL diagnostics for Rpath models
#'
#' Implements standard PREBAL diagnostics following:
#' Link (2010) and Heymans et al. (2016)
#'
#' @param rpath.obj A balanced Rpath model object
#'
#' @return A list containing structural metrics, slope diagnostics,
#'         vital rate checks, and flagged groups
#' @export
run_prebal <- function(rpath.obj) {
  
  if (!requireNamespace("Rpath", quietly = TRUE)) {
    stop("Package 'Rpath' required.")
  }
  
  # Extract Ecopath-style output
  model <- as.data.frame(Rpath::write.Rpath(rpath.obj))
  
  # Living groups only (type < 2 in Rpath convention)
  living <- model[model$type < 2, ]
  

  # 1. BIOMASS STRUCTURE ----------
  
  biomass_span <- log10(max(living$Biomass, na.rm = TRUE)) -
    log10(min(living$Biomass[living$Biomass > 0], na.rm = TRUE))
  
  biomass_lm <- lm(log10(Biomass) ~ TL, data = living)
  biomass_slope <- as.numeric(coef(biomass_lm)[2])
  biomass_se    <- as.numeric(summary(biomass_lm)$coef[2, 2])
  

  # 2. TROPHIC AGGREGATED BIOMASS SLOPE --------------
  # (Biomass summed within TL bins)

  
  TL_breaks <- seq(floor(min(living$TL)),
                   ceiling(max(living$TL)),
                   by = 0.5)
  
  living$TL_bin <- cut(living$TL,
                       breaks = TL_breaks,
                       include.lowest = TRUE)
  
  TL_sum <- aggregate(Biomass ~ TL_bin, data = living, sum)
  TL_mid <- sapply(strsplit(as.character(TL_sum$TL_bin), ","), function(x) {
    nums <- gsub("\\[|\\]|\\(|\\)", "", x)
    vals <- as.numeric(strsplit(nums, ",")[[1]])
    mean(vals)
  })
  
  TL_lm <- lm(log10(TL_sum$Biomass) ~ TL_mid)
  TL_slope <- coef(TL_lm)[2]
  TL_se    <- summary(TL_lm)$coef[2, 2]
  

  # 3. VITAL RATE SLOPES -------------------
  
  # P/B slope
  pb_lm <- lm(log10(PB) ~ TL, data = living)
  pb_slope <- coef(pb_lm)[2]
  pb_se    <- summary(pb_lm)$coef[2, 2]
  
  # Q/B slope (exclude primary producers)
  qb_data <- living[living$QB > 0 & living$Group != "Phytoplankton", ]
  qb_lm <- lm(log10(QB) ~ TL, data = qb_data)
  qb_slope <- coef(qb_lm)[2]
  qb_se    <- summary(qb_lm)$coef[2, 2]
  

  # 4. GROSS EFFICIENCY CHECK (P/Q) ---------------------
  # Expected ~0.1–0.3
  
  consumer_data <- living[living$QB > 0 & living$type == 0, ]
  GE_flag <- consumer_data[consumer_data$GE < 0.1 | consumer_data$GE > 0.3, ]
  
  # 5. ECOLOGICAL EFFICIENCY CHECK (EE <= 1) --------------------

  
  EE_flag <- living[living$EE > 1,
                    c("Group", "EE")]
  

  # 6. FISHING MORTALITY CHECK -------------------
  # F should not exceed total mortality (P/B)
  
  living$Fmort <- ifelse(living$Biomass > 0,
                         living$Removals / living$Biomass,
                         NA)
  
  F_flag <- living[living$Fmort > living$PB,
                   c("Group", "Biomass", "Removals", "Fmort", "PB")]
  
  # 7. SUMMARY INTERPRETATION FLAGS ------------
  
  structure_flags <- list(
    biomass_span_ok = biomass_span >= 5,
    biomass_slope_expected_range = biomass_slope < -0.3 & biomass_slope > -1.5,
    TL_slope_expected_range = TL_slope < -0.1 & TL_slope > -1.5,
    pb_declines_with_TL = pb_slope < 0,
    qb_declines_with_TL = qb_slope < 0
  )
  

  # RETURN CLEAN STRUCTURED OUTPUT ------------------

  
  return(list(
    structure = list(
      biomass_span = biomass_span,
      biomass_slope = biomass_slope,
      biomass_slope_se = biomass_se,
      TL_slope = TL_slope,
      TL_slope_se = TL_se
    ),
    vital_rates = list(
      pb_slope = pb_slope,
      pb_slope_se = pb_se,
      qb_slope = qb_slope,
      qb_slope_se = qb_se
    ),
    flags = list(
      GE_out_of_range = GE_flag,
      EE_gt_1 = EE_flag,
      fishing_exceeds_PB = F_flag
    ),
    interpretation = structure_flags
  ))
}
