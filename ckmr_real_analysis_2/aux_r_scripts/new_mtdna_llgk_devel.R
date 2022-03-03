# =================================================================
# Development of new CIFF lglk. Start in R and then transfer
# Author: Luke Lloyd-Jones 
# Date started: 13/12/2021
# Date updated: 16/12/2021
# =================================================================
# Libraries needed
suppressPackageStartupMessages({
  library(atease)   
  library(mvbutils)
  library(dplyr)
  library(offarray)
  library(reshape)
  library(RColorBrewer)
  library(TMB)
  library(ggplot2)
})

# Source functions required for simulation and analysis
source("functions/ciff_llk_mtdna.R")

#------------------------------------------------------------------------
# Read back in the kin sets and the pop sample and prime the combinations
# and kin for fitting.
# All kin and combinations have been generated in previous setup in the
# ../kin
# Loop over all possible comparisons for HSPs and just ask are they HSP
# but adding the PHSP and MHSP probabilities together.
# Then we have an additional loop over over those pairs for which the 
# additional data are available i.e., not for the missing ones which looks at 
# the conditional lglk of the pairwise datum given that the pair is kin of 
# whatever specific type (and the pair's other covariates).
# ------------------------------------------------------------------------


# --------------------------------------
# OK now with the special mtDNA loop
# --------------------------------------
mop.kin <- get(load("../kin_and_combos/comp_kin_out/mop_kin.Rdata"))
fop.kin <- get(load("../kin_and_combos/comp_kin_out/fop_kin.Rdata"))
sib.kin <- get(load("../kin_and_combos/comp_kin_out/sibling_kin.Rdata"))

Yad_range <- as.numeric(colnames(mop.kin$mop.comp))
Bju_range <- as.numeric(colnames(sib.kin$sib.comp))

pars <- c()
# With just deffo sibling probabilities
fit2 <- nlminb(pars, 
               ciff_ckmr_loglik_mtdna, 
               lower = c(1,     -20,     0.01, -10), 
               upper = c(20000, 20, 1,     10))

# Mitochondrial loop
defo.sibs <- read.csv("../mtDNA/results/all_siblings_plod34_mtdna_v_haplo.csv")
# Add a HSP and FSP index onto this as we need this in the clean up loop.
defo.sibs$FSP <- 0
defo.sibs$FSP[defo.sibs$PLOD > 100] <- 1
# Need to exclude those that have no haplo info
defo.sibs <- defo.sibs[-which(is.na(defo.sibs$HAP1) | 
                              is.na(defo.sibs$HAP2)), ]

#diag(sib.kin$sib.comp) <- 0
#diag(sib.kin$sib.kin)  <- 0
defo.sibs <- defo.sibs[-which(defo.sibs$Ch_Gap == 0), ]

data      <- list(comb_mop  = as.array(mop.kin$mop.comp),    
                  kin_mop   = as.array(mop.kin$mop.kin),
                  comb_fop  = as.array(fop.kin$ppop.comp),    
                  kin_fop   = as.array(fop.kin$ppop.kin),
                  comb_sib  = as.array(sib.kin$sib.comp),       
                  kin_sib   = as.array(sib.kin$sib.kin),
                  symin     = min(Yad_range), 
                  symax     = max(Yad_range), 
                  bymin     = min(Bju_range),  
                  bymax     = max(Bju_range), 
                  lethal    = 0, 
                 surv_prior = 1,
               fudge_factor = 1  - 0.081,
                      Mu_Z1 = 0.565, 
                      SD_Z1 = 0.07,
                  sib_pairs = defo.sibs)

pars <- c(log(500), 
          log(0.8 / (1 - 0.8)), 
          0.5,
          0.0)
ciff_ckmr_loglik_mtdna(pars)
pars <- c(log(500), 
          log(0.8 / (1 - 0.8)), 
          0.5,
          0.0,
          0.1)

pars <- c(log(1000), 
          log(0.8 / (1 - 0.8)), 
          0.62,
          0.0,
          0.5)

pars <- c(log(1000), 
          log(0.8 / (1 - 0.8)), 
          0.62,
          0.0,
          0.5)

#pars <- c(7.14474129,  0.44636374, 0.6, -0.09325861,  8.11912435)
fit2 <- nlminb(pars, 
               NEG(ciff_ckmr_loglik_mtdna), 
               lower = c(1,  -20,  0, -10, -20), 
               upper = c(20,  20,  1,  10,  20))
fit2

pars <- fit2$par


exp.mats <- ciff_ckmr_loglik_mtdna(pars, ret_mats = 1)[[2]]


# -----------------------------
# Matrix diagonal for reporting
# -----------------------------
# Maternal
get.ods <- function(A, ofst)
{
  A[row(A) == (col(A) - ofst)]
}
sum(get.ods(mop.kin$mop.comp, -1)) + sum(get.ods(mop.kin$mop.comp, -3)) # 1 year diff
sum(get.ods(mop.kin$mop.comp, 0)) + sum(get.ods(mop.kin$mop.comp, -4)) # 2 year diff
sum(get.ods(mop.kin$mop.comp, 1)) + sum(get.ods(mop.kin$mop.comp, -5)) # 3 year diff
sum(get.ods(mop.kin$mop.comp, 2))
sum(get.ods(mop.kin$mop.comp, 3))
sum(get.ods(mop.kin$mop.comp, 4))

sum(get.ods(mop.kin$mop.kin, -1)) + sum(get.ods(mop.kin$mop.kin, -3)) # 1 year diff
sum(get.ods(mop.kin$mop.kin, 0)) + sum(get.ods(mop.kin$mop.kin, -4)) # 2 year diff
sum(get.ods(mop.kin$mop.kin, 1)) + sum(get.ods(mop.kin$mop.kin, -5)) # 3 year diff
sum(get.ods(mop.kin$mop.kin, 2)) # 4
sum(get.ods(mop.kin$mop.kin, 3)) # 5
sum(get.ods(mop.kin$mop.kin, 4)) # 6

# Paternal

# Combinations
sum(get.ods(fop.kin$ppop.comp, -1)) + sum(get.ods(fop.kin$ppop.comp, -3)) # 1 year diff
sum(get.ods(fop.kin$ppop.comp, 0))  + sum(get.ods(fop.kin$ppop.comp, -4)) # 2 year diff
sum(get.ods(fop.kin$ppop.comp, 1))  + sum(get.ods(fop.kin$ppop.comp, -5)) # 3 year diff
sum(get.ods(fop.kin$ppop.comp, 2)) # 4
sum(get.ods(fop.kin$ppop.comp, 3)) # 5
sum(get.ods(fop.kin$ppop.comp, 4)) # 6
sum(fop.kin$ppop.comp)
# Kin
sum(get.ods(fop.kin$ppop.kin, -1)) + sum(get.ods(fop.kin$ppop.kin, -3)) # 1 year diff
sum(get.ods(fop.kin$ppop.kin, 0))  + sum(get.ods(fop.kin$ppop.kin, -4)) # 2 year diff
sum(get.ods(fop.kin$ppop.kin, 1))  + sum(get.ods(fop.kin$ppop.kin, -5)) # 3 year diff
sum(get.ods(fop.kin$ppop.kin, 2)) # 4
sum(get.ods(fop.kin$ppop.kin, 3)) # 5
sum(get.ods(fop.kin$ppop.kin, 4)) # 6
sum(fop.kin$ppop.kin)

# --------
# Siblings
# --------

# Combinations
sum(get.ods(sib.kin$sib.comp, 0))  
sum(get.ods(sib.kin$sib.comp, 1)) 
sum(get.ods(sib.kin$sib.comp, 2))
sum(get.ods(sib.kin$sib.comp, 3))
sum(get.ods(sib.kin$sib.comp, 4))
sum(get.ods(sib.kin$sib.comp, 5))
sum(sib.kin$sib.comp)
# Kin 

sum(get.ods(sib.kin$sib.kin, 0))  
sum(get.ods(sib.kin$sib.kin, 1)) 
sum(get.ods(sib.kin$sib.kin, 2)) # 4
sum(get.ods(sib.kin$sib.kin, 3)) # 5
sum(get.ods(sib.kin$sib.kin, 4)) # 6
sum(get.ods(sib.kin$sib.kin, 5)) # 6
sum(sib.kin$sib.kin)



# ------------------------------
# Expected matrices
# ------------------------------

# MOPS
signif(sum(get.ods(exp.mats$mops, -1)) + sum(get.ods(exp.mats$mops, -3)), 3) # 1 year diff
signif(sum(get.ods(exp.mats$mops, 0))  + sum(get.ods(exp.mats$mops, -4)), 3) # 2 year diff
signif(sum(get.ods(exp.mats$mops, 1))  + sum(get.ods(exp.mats$mops, -5)), 3) # 3 year diff
signif(sum(get.ods(exp.mats$mops, 2)), 3) # 4
signif(sum(get.ods(exp.mats$mops, 3)), 3) # 5
signif(sum(get.ods(exp.mats$mops, 4)), 3) # 6
sum(exp.mats$mops)

# FOPS
signif(sum(get.ods(exp.mats$fops, -1)) + sum(get.ods(exp.mats$fops, -3)), 3) # 1 year diff
signif(sum(get.ods(exp.mats$fops, 0))  + sum(get.ods(exp.mats$fops, -4)), 3) # 2 year diff
signif(sum(get.ods(exp.mats$fops, 1))  + sum(get.ods(exp.mats$fops, -5)), 3) # 3 year diff
signif(sum(get.ods(exp.mats$fops, 2)), 3) # 4
signif(sum(get.ods(exp.mats$fops, 3)), 3) # 5
signif(sum(get.ods(exp.mats$fops, 4)), 3) # 6
sum(exp.mats$fops)


# Full siblings
signif(sum(get.ods(exp.mats$fsps, 1)), 3)
signif(sum(get.ods(exp.mats$fsps, 2)), 3) # 2
signif(sum(get.ods(exp.mats$fsps, 3)), 3) # 3
signif(sum(get.ods(exp.mats$fsps, 4)), 3) # 4
signif(sum(get.ods(exp.mats$fsps, 5)), 3) # 5
sum(exp.mats$fsps)

# maternal siblings
signif(sum(get.ods(exp.mats$mhsps, 1)), 3)
signif(sum(get.ods(exp.mats$mhsps, 2)), 3) # 2
signif(sum(get.ods(exp.mats$mhsps, 3)), 3) # 3
signif(sum(get.ods(exp.mats$mhsps, 4)), 3) # 4
signif(sum(get.ods(exp.mats$mhsps, 5)), 3) # 5
sum(exp.mats$mhsps)

13+22+4
# paternal siblings
signif(sum(get.ods(exp.mats$fhsps , 1)), 3)
signif(sum(get.ods(exp.mats$fhsps, 2)), 3) # 2
signif(sum(get.ods(exp.mats$fhsps, 3)), 3) # 3
signif(sum(get.ods(exp.mats$fhsps, 4)), 3) # 4
signif(sum(get.ods(exp.mats$fhsps, 5)), 3) # 5
sum(exp.mats$fhsps)


# Table the observed data

defo.sibs.fsps <- defo.sibs[which(defo.sibs$FSP == 1), ]

defo.sibs.mhsps <- defo.sibs[which(defo.sibs$SAME_HAP == 1 &
                                  defo.sibs$FSP != 1), ]

defo.sibs.fhsps <- defo.sibs[which(defo.sibs$SAME_HAP != 1 &
                                   defo.sibs$FSP != 1), ]

obs <- c(0, 0, 1, 0, 0, 0)
exp <- c(1.89, 1.23, 0.554, 0.256, 0.029, 0.018)
x <- chisq.test(obs, exp)


(0 - 1.89)^2/1.89 +  (0-1.23)^2/1.23    +    (1-0.554)^2/0.554 + (0-0.256)^2/0.256  +  (0 - 0.029)^2/0.029 + (0 - 0.018)^2/0.018
pchisq(3.782054, df = 5)

obs <- c(3,  2,   0,       0,  0, 0)
exp <- c(8.11, 3.8, 1.47, 0.635, 0.066, 0.044)
scr <- sum((obs-exp)^2/exp)
1-pchisq(scr, 5)


obs <- c(6 ,    6,   3,  2, 0)
exp <-  c(8.08,  4.06,   1.31   , 0.294, 0.0655)
scr <- sum((obs-exp)^2/exp)
1-pchisq(scr, 4)



obs <- c(14 ,    9,   2   , 1  , 0)
exp <- c(14.2,  6.39,   1.86   , 0.374, 0.076)
scr <- sum((obs-exp)^2/exp)
1-pchisq(scr, 4)

14/14.2/20718        &  9/6.39/14660    &   2/1.86/6864    & 1/0.374/2250     & 0/0.076/759



