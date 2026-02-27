#' Plot ranked biomasses for balanced Rpath model
#'
#' Produces horizontal bar graph of living-group biomasses ranked
#' from largest (top) to smallest (bottom).
#'
#' @param rpath.obj A balanced Rpath model object
#' @param log_scale Logical; if TRUE uses log10 biomass
#' @param save Logical; if TRUE saves plot to file
#' @param filename Output filename if save = TRUE
#'
#' @return Invisibly returns ranked biomass data.frame
#' @export

plot_balanced_biomass <- function(rpath.obj,
                                  log_scale = TRUE,
                                  save = FALSE,
                                  filename = "balanced_biomass_ranked.png") {
  
  if (!requireNamespace("Rpath", quietly = TRUE)) {
    stop("Package 'Rpath' required.")
  }
  
  biomass_df <- diagnostics_balanced_biomass(rpath.obj)
  
  if (log_scale) {
    biomass_values <- log10(biomass_df$Biomass)
    xlab <- "log10 Biomass"
  } else {
    biomass_values <- biomass_df$Biomass
    xlab <- "Biomass"
  }
  
  # Reverse order so largest appears at top
  biomass_df <- biomass_df[rev(seq_len(nrow(biomass_df))), ]
  biomass_values <- rev(biomass_values)
  
  # Color by trophic level
  tl_scaled <- (biomass_df$TL - min(biomass_df$TL)) /
    (max(biomass_df$TL) - min(biomass_df$TL))
  
  bar_cols <- colorRampPalette(c("steelblue", "darkred"))(100)[
    as.numeric(cut(tl_scaled, breaks = 100))
  ]
  
  n_groups <- nrow(biomass_df)
  
  if (save) {
    png(filename,
        width = 1200,
        height = max(1600, n_groups * 30),
        res = 150)
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar = c(4, 12, 3, 2))
  
  barplot(height = biomass_values,
          names.arg = biomass_df$Group,
          horiz = TRUE,
          las = 1,
          col = bar_cols,
          border = NA,
          xlab = xlab,
          main = "Ranked Living-Group Biomass",
          cex.names = max(0.5, 1 - (n_groups / 100)))
  
  if (save) dev.off()
  
  invisible(biomass_df)
}
