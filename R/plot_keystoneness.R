# R/plot_keystoneness.R

#' Plot keystoneness index values
#'
#' Creates a barplot of keystoneness index values (Libralato et al. 2006)
#' for functional groups in an Rpath model. Groups are ordered by KS value,
#' and positive values (potential keystone groups) are highlighted.
#'
#' @param ks_results A data frame returned by [keystoneness()],
#'   with at least columns: `group` and `keystoneness`.
#' @param top_n Optional. Highlight the top N groups with the highest KS values.
#'   Default = 5.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#'   ks <- neusRpath::keystoneness("GB")
#'   neusRpath::plot_keystoneness(ks)
#' }
plot_keystoneness <- function(ks_results, top_n = 5) {
  stopifnot(all(c("group", "keystoneness") %in% names(ks_results)))
  
  library(ggplot2)
  library(dplyr)
  
  # Order groups by KS
  plot_df <- ks_results %>%
    arrange(keystoneness) %>%
    mutate(group = factor(group, levels = group),
           highlight = ifelse(rank(-keystoneness) <= top_n, "Top", "Other"))
  
  # Plot
  p <- ggplot(plot_df, aes(x = group, y = keystoneness, fill = highlight)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Top" = "firebrick", "Other" = "grey70")) +
    labs(x = "Functional Group", y = "Keystoneness (Libralato et al. 2006)",
         fill = paste("Top", top_n)) +
    theme_minimal(base_size = 12)
  
  return(p)
}
