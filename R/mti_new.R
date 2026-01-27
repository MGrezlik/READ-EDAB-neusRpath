#' Mixed Trophic Impact (Ecopath-consistent)
#'
#' Computes the Mixed Trophic Impact (MTI) matrix for an Rpath model,
#' reproducing Ecopath behavior including undocumented legacy steps.
#'
#' @param Rpath A balanced Rpath model object
#' @param Rpath.params The corresponding unbalanced/balanced parameter object
#'
#' @return A square MTI matrix with rows and columns equal to all groups
#' @export
MTI_new <- function(Rpath, Rpath.params) {
  
  ## ------------------------------------------------------------------
  ## Basic dimensions and indices
  ## ------------------------------------------------------------------
  
  NUM_LIVING <- Rpath$NUM_LIVING
  NUM_DEAD   <- Rpath$NUM_DEAD
  NUM_GROUPS <- Rpath$NUM_GROUPS
  N_BIO      <- NUM_LIVING + NUM_DEAD
  
  group_names <- Rpath$Group
  
  ## ------------------------------------------------------------------
  ## Initialize matrices
  ## ------------------------------------------------------------------
  
  AA  <- matrix(0, NUM_GROUPS, NUM_GROUPS,
                dimnames = list(group_names, group_names))
  
  DCC <- matrix(0, NUM_GROUPS, NUM_GROUPS,
                dimnames = list(group_names, group_names))
  
  ## ------------------------------------------------------------------
  ## Diet composition matrix (DCC)
  ##   - remove Group column explicitly
  ##   - remove Import row
  ##   - add detrital columns
  ## ------------------------------------------------------------------
  
  diet_df <- Rpath.params$diet
  
  if (!("Group" %in% names(diet_df))) {
    stop("Rpath.params$diet must contain a 'Group' column")
  }
  
  diet_mat <- as.matrix(
    diet_df[, setdiff(names(diet_df), "Group")]
  )
  
  # Remove Import row (assumed last row, Ecopath convention)
  diet_mat <- diet_mat[1:(nrow(diet_mat) - 1), , drop = FALSE]
  
  # Fill biological portion of DCC
  DCC[1:N_BIO, 1:NUM_LIVING] <- diet_mat
  
  # Add detritus columns (zeros)
  if (NUM_DEAD > 0) {
    DCC[1:N_BIO, (NUM_LIVING + 1):N_BIO] <- 0
  }
  
  ## ------------------------------------------------------------------
  ## Fishing columns in DCC (landings only; Ecopath behavior)
  ## ------------------------------------------------------------------
  
  if (!is.null(Rpath$Landings)) {
    DCC[, (N_BIO + 1):NUM_GROUPS] <- Rpath$Landings
  }
  
  ## ------------------------------------------------------------------
  ## Normalize DCC by column sums
  ## ------------------------------------------------------------------
  
  DCC[is.na(DCC)] <- 0
  col_sums <- colSums(DCC)
  
  DCC <- sweep(DCC, 2, col_sums, "/")
  DCC[!is.finite(DCC)] <- 0
  
  ## ------------------------------------------------------------------
  ## Flow matrix AA (consumption + fishing + discards)
  ## ------------------------------------------------------------------
  
  for (i in seq_len(N_BIO)) {
    for (j in seq_len(NUM_LIVING)) {
      AA[i, j] <- Rpath$Biomass[j] * Rpath$QB[j] * DCC[i, j]
    }
  }
  
  # Fishing removals (landings + discards)
  if (!is.null(Rpath$Landings)) {
    AA[, (N_BIO + 1):NUM_GROUPS] <- Rpath$Landings
  }
  
  if (!is.null(Rpath$Discards)) {
    AA[, (N_BIO + 1):NUM_GROUPS] <-
      AA[, (N_BIO + 1):NUM_GROUPS] + Rpath$Discards
  }
  
  # Discards routed to detritus
  if (!is.null(Rpath$Discards) && !is.null(Rpath$DetFate)) {
    fishing_det <- colSums(Rpath$Discards) *
      Rpath$DetFate[(N_BIO + 1):NUM_GROUPS, , drop = FALSE]
    
    AA[(N_BIO + 1):NUM_GROUPS,
       (NUM_LIVING + 1):N_BIO] <- fishing_det
  }
  
  ## ------------------------------------------------------------------
  ## Detritus is not a predator
  ## ------------------------------------------------------------------
  
  AA[, (NUM_LIVING + 1):N_BIO] <- 0
  
  ## ------------------------------------------------------------------
  ## Fraction of consumption (FC)
  ## ------------------------------------------------------------------
  
  row_sums <- rowSums(AA)
  FC <- sweep(AA, 1, row_sums, "/")
  FC[!is.finite(FC)] <- 0
  
  FC_T <- t(FC)
  
  ## ------------------------------------------------------------------
  ## Direct impact matrix
  ## ------------------------------------------------------------------
  
  MTI_direct <- DCC - FC_T
  
  ## ------------------------------------------------------------------
  ## Total MTI: (I − MTI)^−1 − I
  ## ------------------------------------------------------------------
  
  MTI_diag <- diag(NUM_GROUPS) - MTI_direct
  MTI_inv  <- MASS::ginv(MTI_diag)
  
  diag(MTI_inv) <- diag(MTI_inv) - 1
  
  ## ------------------------------------------------------------------
  ## Zero detritus self-impacts (Ecopath legacy behavior)
  ## ------------------------------------------------------------------
  
  if (NUM_DEAD > 0) {
    det_idx <- (NUM_LIVING + 1):N_BIO
    MTI_inv[cbind(det_idx, det_idx)] <- 0
  }
  
  ## ------------------------------------------------------------------
  ## Final formatting
  ## ------------------------------------------------------------------
  
  colnames(MTI_inv) <- group_names
  rownames(MTI_inv) <- group_names
  
  return(MTI_inv)
}


#############################################################
# TESTS

# load(here::here('data/GB.rda'))
# load(here::here('data/GB_balanced_params.rda'))
# 
# MTI_out <- MTI_new(GB, GB_balanced_params)
# 
# MTI_Rpath <- Rpath::MTI(GB, GB_balanced_params)
# 
# final_diff <- MTI_Rpath - MTI_out
# 
# max(final_diff)
# 
# # From Valls et al
# epsilon_mat <- MTI_out[1:GB$NUM_LIVING, 1:GB$NUM_LIVING]
# diag(epsilon_mat) <- 0
# epsilon_i <- sqrt(rowSums(epsilon_mat * epsilon_mat))
# 
# p_i <- GB$Biomass[1:GB$NUM_LIVING]/sum(GB$Biomass[1:GB$NUM_LIVING])
# 
# KS_1 <- log10(epsilon_i * (1-p_i))
# KS_2 <- log10(epsilon_i * (1/p_i))
# # Biomass Rank in descending order
# drank <- GB$NUM_LIVING - rank(GB$Biomass[1:GB$NUM_LIVING]) + 1
# KS_3 <- log10(epsilon_i * drank)
# 
# Rel.Tot <- (rowSums((MTI_out[1:GB$NUM_LIVING, 1:GB$NUM_LIVING])))  
# 
# #test <- cbind(ewe_ks,KS_1,KS_2,Rel.Tot,KS_3)
# 
# KS_out <- cbind(KS_1,KS_2,KS_3)
# 
# write.csv(KS_out_GB, "GB_keystones.csv", row.names=T)
# write.csv(MTI_out_GB, "GB_MTI_matrix.csv", row.names=T)
# 
# 
# t1 <- abs(epsilon_mat)
# t2 <- t(t(t1)/rowSums(t(t1)))
# 
# rowSums(t2)