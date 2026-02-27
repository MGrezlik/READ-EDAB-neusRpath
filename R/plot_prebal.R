#' Plot PREBAL diagnostics for Rpath models
#'
#' Produces standard PREBAL diagnostic plots:
#' 1) Biomass vs TL
#' 2) Aggregated biomass by TL bins
#' 3) P/B vs TL
#' 4) Q/B vs TL
#'
#' @param rpath.obj A balanced Rpath model object
#' @param save Logical. If TRUE, saves plot to file.
#' @param filename Output filename if save = TRUE
#'
#' @return Invisibly returns list of fitted model objects
#' @export

plot_prebal <- function(rpath.obj,
                        save = FALSE,
                        filename = "prebal_diagnostics.png") {
  
  if (!requireNamespace("Rpath", quietly = TRUE)) {
    stop("Package 'Rpath' required.")
  }
  
  model <- as.data.frame(Rpath::write.Rpath(rpath.obj))
  living <- model[model$type < 2 & model$Biomass > 0, ]
  
  # TL bins for aggregated biomass
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
  
  # Fit models
  biomass_lm <- lm(log10(Biomass) ~ TL, data = living)
  TL_lm      <- lm(log10(TL_sum$Biomass) ~ TL_mid)
  pb_lm      <- lm(log10(PB) ~ TL, data = living)
  qb_data    <- living[living$QB > 0 & living$Group != "Phytoplankton", ]
  qb_lm      <- lm(log10(QB) ~ TL, data = qb_data)
  
  if (save) {
    png(filename, width = 1400, height = 1200, res = 150)
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(2,2), mar = c(4,4,3,1))
  
  # 1. Biomass vs TL -------------
  plot(living$TL,
       log10(living$Biomass),
       pch = 16,
       xlab = "Trophic Level",
       ylab = "log10 Biomass",
       main = "Biomass vs Trophic Level")
  
  abline(biomass_lm, lwd = 2)
  
  mtext(paste0("Slope = ",
               round(unname(coef(biomass_lm)[2]),3),
               ", R² = ",
               round(summary(biomass_lm)$r.squared,2)),
        side = 3, line = -1.5, adj = 1, cex = 0.8)
  
  # 2. Aggregated Biomass --------------
  plot(TL_mid,
       log10(TL_sum$Biomass),
       pch = 16,
       xlab = "Trophic Level (bin midpoints)",
       ylab = "log10 Summed Biomass",
       main = "Aggregated Biomass by TL")
  
  abline(TL_lm, lwd = 2)
  
  mtext(paste0("Slope = ",
               round(unname(coef(TL_lm)[2]),3),
               ", R² = ",
               round(summary(TL_lm)$r.squared,2)),
        side = 3, line = -1.5, adj = 1, cex = 0.8)
  
  # 3. P/B vs TL --------------
  plot(living$TL,
       log10(living$PB),
       pch = 16,
       xlab = "Trophic Level",
       ylab = "log10 P/B",
       main = "Production/Biomass vs TL")
  
  abline(pb_lm, lwd = 2)
  
  mtext(paste0("Slope = ",
               round(unname(coef(pb_lm)[2]),3),
               ", R² = ",
               round(summary(pb_lm)$r.squared,2)),
        side = 3, line = -1.5, adj = 1, cex = 0.8)
  
  # 4. Q/B vs TL -----------------
  plot(qb_data$TL,
       log10(qb_data$QB),
       pch = 16,
       xlab = "Trophic Level",
       ylab = "log10 Q/B",
       main = "Consumption/Biomass vs TL")
  
  abline(qb_lm, lwd = 2)
  
  mtext(paste0("Slope = ",
               round(unname(coef(qb_lm)[2]),3),
               ", R² = ",
               round(summary(qb_lm)$r.squared,2)),
        side = 3, line = -1.5, adj = 1, cex = 0.8)
  
  if (save) dev.off()
  
  invisible(list(
    biomass_lm = biomass_lm,
    TL_lm = TL_lm,
    pb_lm = pb_lm,
    qb_lm = qb_lm
  ))
}