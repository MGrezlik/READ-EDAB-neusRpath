# MTG 01/27/2026
# keystoneness values calculated as described in the literature do not align with values
# calculated in Ecopath.
# Kerim found old code from Villy which corrects the MTI function
# MTI_new function is copied from his rpath_notes repo
# That repo is public for now but will need to be private

library(Rpath)

MTI_new <- function(Rpath, Rpath.params){
  
  # Useful indices
  N_BIO <- Rpath$NUM_LIVING+Rpath$NUM_DEAD
  
  # Make empty matrices
  AA <- matrix(0, Rpath$NUM_GROUPS, Rpath$NUM_GROUPS)
  colnames(AA) <- Rpath$Group
  rownames(AA) <- Rpath$Group
  DCC <- AA
  
  # Diet matrix removing import row, group column, adding detrital columns
  DCC[1:N_BIO,1:N_BIO] <- cbind(as.matrix(Rpath.params$diet[,-"Group"])[1:(dim(Rpath.params$diet)[1]-1),],
                                matrix(0, nrow=dim(Rpath.params$diet)[1]-1, ncol=Rpath$NUM_DEAD))
  
  # Landings only not discards in DCC (discards never included at all in this??)
  DCC[,(N_BIO+1):Rpath$NUM_GROUPS] <- Rpath$Landings 
  
  # Normalizing and cleaning DCC (diet composition and landings)
  DCC[is.na(DCC)] <- 0
  DCC <- t(t(DCC)/colSums(DCC))
  DCC[is.nan(DCC) | is.infinite(DCC)] <- 0
  
  # Living total flows (and flows FROM detritus via diets, but not TO detritus)
  for (i in 1:N_BIO){
    for (j in 1:Rpath$NUM_LIVING){
      AA[i,j] <- Rpath$Biomass[j] * Rpath$QB[j] * DCC[i,j]
    }
  }
  
  # These lines not needed because flow into detritus in AA matrix is 0'd out below?
  # Keeping in comments for reference.
  #  detflow <- (Rpath$PB*Rpath$Biomass*(1-Rpath$EE)) + Rpath$QB*Rpath$Biomass*Rpath$Unassim
  #  detmat <- detflow * Rpath$DetFate
  #  AA[,(Rpath$NUM_LIVING+1):N_BIO] <- detmat
  # Should Discards be here?  Seems to work but check later?  
  AA[,(N_BIO+1):Rpath$NUM_GROUPS] <- Rpath$Landings + Rpath$Discards 
  fishing_det <- colSums(Rpath$Discards) * Rpath$DetFate[(N_BIO+1):Rpath$NUM_GROUPS,]
  AA[(N_BIO+1):Rpath$NUM_GROUPS, (Rpath$NUM_LIVING+1):N_BIO] <- fishing_det
  
  # "Detritus is not considered a predator" - remove flows to detritus before normalizing
  AA[,(Rpath$NUM_LIVING+1):N_BIO] <- 0
  
  # Normalize FC and and then flip so that predators are in columns    
  FC <- AA/rowSums(AA)
  FC[is.nan(FC) | is.infinite(FC)] <- 0
  FC_T <- t(FC)
  
  # Direct impacts matrix
  MTI <- (DCC - FC_T)
  
  # From Ulanowicz paper, find inverse of ([Identity] - MTI)
  MTI_diag <- diag(Rpath$NUM_GROUPS) - MTI
  
  # Invert the matrix    
  MTI_inv <- MASS::ginv((MTI_diag))
  
  # Now subtract 1 off the diagonal matrix  
  diag(MTI_inv) <- diag(MTI_inv) - 1
  
  # The final(?) piece, hidden in Villy's '03 code - 0 out the detritial Identity 
  for (i in (Rpath$NUM_LIVING+1):N_BIO){
    MTI_inv[i,i] <- 0
  }
  
  colnames(MTI_inv) <- Rpath$Group
  rownames(MTI_inv) <- Rpath$Group
  
  return(MTI_inv)  
}

#############################################################
# TESTS

load(here::here('data'))

MTI_out <- MTI_new(bal, unbal)

inmat <- read.csv(mat_file)  
rownames(inmat) <- inmat$GROUP
imat <- as.matrix(inmat[,-1])
final_diff <- imat - MTI_out

max(final_diff)

##############################################################
#
#this_model <- "xml_examples/WGOA_9May2025.eiixml"
#mat_file <- "WGOA_mat.csv"  
#ks_file  <- "WGOA_ks.csv"  
#unbal   <- rpath.stanzas(create.rpath.from.eiixml(this_model))

unbal <- rpath.stanzas(read.rpath.params(
  modfile   = "rpath_examples/ebs_full_base_oct2024_fisheries_area_fix.csv",
  dietfile  = "rpath_examples/ebs_full_diet_apr2025_alt_invert_diets.csv",
  pedfile   = NA,
  stanzagroupfile = "rpath_examples/EBS_full_stanza_groups_feb2023.csv",
  stanzafile      = "rpath_examples/EBS_full_stanzas_mar2023.csv"
))


unbal <- rpath.stanzas(read.rpath.params(
  modfile   = "rpath_examples/GOA_full_base.csv",
  dietfile  = "rpath_examples/GOA_full_diet.csv",
  pedfile   = NA,
  stanzagroupfile = "rpath_examples/GOA_full_stanza_groups.csv",
  stanzafile      = "rpath_examples/GOA_full_stanzas.csv"
))  

bal     <- rpath(unbal)  
MTI_out <- MTI_new(bal, unbal)  

# Comparison file
#ewe_ks <- read.csv(ks_file)

# From Valls et al
epsilon_mat <- MTI_out[1:bal$NUM_LIVING, 1:bal$NUM_LIVING]
diag(epsilon_mat) <- 0
epsilon_i <- sqrt(rowSums(epsilon_mat * epsilon_mat))

p_i <- bal$Biomass[1:bal$NUM_LIVING]/sum(bal$Biomass[1:bal$NUM_LIVING])

KS_1 <- log10(epsilon_i * (1-p_i))
KS_2 <- log10(epsilon_i * (1/p_i))
# Biomass Rank in descending order
drank <- bal$NUM_LIVING - rank(bal$Biomass[1:bal$NUM_LIVING]) + 1
KS_3 <- log10(epsilon_i * drank)

Rel.Tot <- (rowSums((MTI_out[1:bal$NUM_LIVING, 1:bal$NUM_LIVING])))  

#test <- cbind(ewe_ks,KS_1,KS_2,Rel.Tot,KS_3)

KS_out <- cbind(KS_1,KS_2,KS_3)

write.csv(KS_out, "GOA_keystones.csv", row.names=T)
write.csv(MTI_out, "GOA_MTI_matrix.csv", row.names=T)


t1 <- abs(epsilon_mat)
t2 <- t(t(t1)/rowSums(t(t1)))

rowSums(t2)