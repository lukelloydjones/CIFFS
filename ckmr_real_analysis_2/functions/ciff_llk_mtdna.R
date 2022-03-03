# ==============================================================
# Try and write CKMR function in R so we can easily diagnose
# problemos with the incorporation of full sibling pairs
# Author: Luke Lloyd-Jones
# Date started: 13/12/2021
# Date updated: 13/12/2021
# ==============================================================


ciff_ckmr_loglik <- function(pars)
{
  N0_f  <- exp(pars[1])
  N0_m  <- exp(pars[1]) * (exp(pars[2]) / (1 + exp(pars[2]))) 
  Z.1   <- pars[3]
  #Z.1   <- 0.6
  R.1   <- pars[4]

  yrs_tot  = data$symax - data$bymin + 1;
  syrs_tot = data$symax - data$symin + 1;
  byrs_tot = data$bymax - data$bymin + 1;
  start_yrs_diff = data$symin - data$bymin;
  y0 <- data$bymin
  # ============================================
  # Female dynamics, survival, and probabilities
  # ============================================
  ## Population dynamics
  years <- seq(data$bymin, data$symax)
  Nf <- offarray( 0, first = data$bymin, 
                     last  = data$symax)
  Nf <- autoloop(yrs = years, {
    N0_f * exp( R.1 * (yrs - y0))
  })
  
  
  # offarray::autoloop() evaluates an expression over all combos of the args
  # ... and returns an array indexed by all combos. Like using nested loops.
  # Indices must match n_comp[] and n_POP[], which are set up by boring_data_prep()
  
  Pr_MOP_noage <- autoloop(Bju = Bju_range, Yad = Yad_range, {
    cumul_surv_1 <- exp( -Z.1 * (Yad - Bju))
    cumul_surv_2 <- exp( -Z.1 * (Bju - Yad))
    (Bju < Yad) * cumul_surv_1 / Nf[ Yad] + (Bju > Yad) * cumul_surv_2 / Nf[ Bju]
    #  Decomposition is:
    #    ( Bju <= Yad ) * # was ad still alive at B[ju] ? Otherwise her ERRO was 0
    #    if she was, then her ERRO is given by:
    #    (N[ Bju] * cumul_surv) / N[ Yad]) /   # Pr[ she was a survivor from the pool of ...
    # ... breeders at Bju, given she was adult at Yad]
    #    N[ Bju]                               # / TRO=Total f adults at B[ju]
    # Then cancel the N[Bju] from top & bottom
  })
  
  if( !all( is.finite( Pr_MOP_noage) & (Pr_MOP_noage >= 0) & (Pr_MOP_noage <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # ============================================
  # Male POP dynamics, survival, and probabilities
  # ============================================
  ## Population dynamics
  
  Nm <- offarray( 0, first = data$bymin, 
                     last  = data$symax)
  Nm <- autoloop(yrs = years, {
    N0_m * exp( R.1 * (yrs - y0))
  })
  
  # # offarray::autoloop() evaluates an expression over all combos of the args
  # # ... and returns an array indexed by all combos. Like using nested loops.
  # # Indices must match n_comp[] and n_POP[], which are set up by boring_data_prep()
  
  Pr_FOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    cumul_surv_1 <- exp( -Z.1 * (Yad - Bju))
    cumul_surv_2 <- exp( -Z.1 * (Bju - Yad))
    (Bju < Yad) * cumul_surv_1 / Nm[ Yad] + (Bju > Yad) * cumul_surv_2 / Nm[ Bju]
    #  Decomposition is:
    #    ( Bju <= Yad ) * # was ad still alive at B[ju] ? Otherwise her ERRO was 0
    #    if she was, then her ERRO is given by:
    #    (N[ Bju] * cumul_surv) / N[ Yad]) /   # Pr[ she was a survivor from the pool of ...
    # ... breeders at Bju, given she was adult at Yad]
    #    N[ Bju]                               # / TRO=Total f adults at B[ju]
    # Then cancel the N[Bju] from top & bottom
  })
  
  if( !all( is.finite( Pr_FOP_noage) & (Pr_FOP_noage >= 0) & (Pr_FOP_noage <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # ============================================
  # Maternal HSPs survival, and probabilities
  # ============================================
  
  Pr_MHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # ORDER is important...
    cumul_surv <- exp( -Z.1 * (B2-B1))
    (B2 >= B1) * cumul_surv / Nf[ B2]
  })
  
  # ============================================
  # Paternal HSPs dynamics, survival, and probabilities
  # ============================================
  # Population dynamics
  
  Pr_PHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # ORDER is important...
    cumul_surv <- exp( -Z.1 * (B2-B1))
    (B2 >= B1) * cumul_surv / Nm[ B2]
  })
  
  if( !all( is.finite( Pr_MHSP) & (Pr_MHSP >= 0) & (Pr_MHSP <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  if( !all( is.finite( Pr_PHSP) & (Pr_PHSP >= 0) & (Pr_PHSP <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # -------------------------
  # Sub to get sibling prob
  # and adjust false-neg rate
  # -------------------------
  mhsp_prob = Pr_MHSP * data$fudge_factor
  fhsp_prob = Pr_PHSP * data$fudge_factor
  Pr_SIB    = mhsp_prob + fhsp_prob
  
  # -----------------
  # POP contributions
  # -----------------
  lglk <- 0
  # MOP
  lglk <- lglk - sum( dbinom(data$kin_mop, size = data$comb_mop, prob = Pr_MOP_noage, log = TRUE))
  # FOP
  lglk <- lglk - sum( dbinom(data$kin_fop, size = data$comb_fop, prob = Pr_FOP_noage, log = TRUE))
  
  # -----------------
  # HSP contributions
  # -----------------
  # SIB
  lglk <- lglk - sum( dbinom(data$kin_sib, size = data$comb_sib, prob = Pr_SIB, log = TRUE))
  
  return(lglk)
}

ciff_ckmr_loglik_mtdna <- function(pars, ret_mats = 0)
{
  # Parameter set up
  N0_f  <- exp(pars[1])
  N0_m  <- exp(pars[1]) * (exp(pars[2]) / (1 + exp(pars[2]))) 
  #Z.1   <- pars[3]
  Z_1   <- pars[3]
  R_1   <- pars[4]
  alpha <- pars[5]
  #R_2   <- pars[6]
  Mu_Z1 <- data$Mu_Z1
  SD_Z1 <- data$SD_Z1
  
  # Data set up
  fneg_rate      <- data$fudge_factor 
  yrs_tot        <- data$symax - data$bymin + 1;
  syrs_tot       <- data$symax - data$symin + 1;
  byrs_tot       <- data$bymax - data$bymin + 1;
  start_yrs_diff <- data$symin - data$bymin;
  y0             <- data$bymin
  sib_pairs      <- data$sib_pairs
  # ============================================
  # Female dynamics, survival, and probabilities
  # ============================================
  ## Population dynamics
  years <- seq(data$bymin, data$symax)
  Nf <- offarray( 0, first = data$bymin, 
                  last  = data$symax)
  Nf <- autoloop(yrs = years, {
    N0_f * exp( R_1 * (yrs - y0))
  })
  
  
  # offarray::autoloop() evaluates an expression over all combos of the args
  # ... and returns an array indexed by all combos. Like using nested loops.
  # Indices must match n_comp[] and n_POP[], which are set up by boring_data_prep()
  
  Pr_MOP_noage <- autoloop(Bju = Bju_range, Yad = Yad_range, {
    cumul_surv_1 <- exp( -Z_1 * (Yad - Bju))
    cumul_surv_2 <- exp( -Z_1 * (Bju - Yad))
    (Bju < Yad) * cumul_surv_1 / Nf[ Yad] + (Bju > Yad) * cumul_surv_2 / Nf[ Bju]
    #  Decomposition is:
    #    ( Bju <= Yad ) * # was ad still alive at B[ju] ? Otherwise her ERRO was 0
    #    if she was, then her ERRO is given by:
    #    (N[ Bju] * cumul_surv) / N[ Yad]) /   # Pr[ she was a survivor from the pool of ...
    # ... breeders at Bju, given she was adult at Yad]
    #    N[ Bju]                               # / TRO=Total f adults at B[ju]
    # Then cancel the N[Bju] from top & bottom
  })
  
  if( !all( is.finite( Pr_MOP_noage) & (Pr_MOP_noage >= 0) & (Pr_MOP_noage <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # ============================================
  # Male POP dynamics, survival, and probabilities
  # ============================================
  ## Population dynamics
  
  Nm <- offarray( 0, first = data$bymin, 
                     last  = data$symax)
  Nm <- autoloop(yrs = years, {
    #N0_m * exp( R_2 * (yrs - y0))
     N0_m * exp( R_1 * (yrs - y0))
  })
  
  # # offarray::autoloop() evaluates an expression over all combos of the args
  # # ... and returns an array indexed by all combos. Like using nested loops.
  # # Indices must match n_comp[] and n_POP[], which are set up by boring_data_prep()
  
  Pr_FOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    cumul_surv_1 <- exp( -Z_1 * (Yad - Bju))
    cumul_surv_2 <- exp( -Z_1 * (Bju - Yad))
    (Bju < Yad) * cumul_surv_1 / Nm[ Yad] + (Bju > Yad) * cumul_surv_2 / Nm[ Bju]
    #  Decomposition is:
    #    ( Bju <= Yad ) * # was ad still alive at B[ju] ? Otherwise her ERRO was 0
    #    if she was, then her ERRO is given by:
    #    (N[ Bju] * cumul_surv) / N[ Yad]) /   # Pr[ she was a survivor from the pool of ...
    # ... breeders at Bju, given she was adult at Yad]
    #    N[ Bju]                               # / TRO=Total f adults at B[ju]
    # Then cancel the N[Bju] from top & bottom
  })
  
  if( !all( is.finite( Pr_FOP_noage) & (Pr_FOP_noage >= 0) & (Pr_FOP_noage <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # ============================================
  # Maternal HSPs survival, and probabilities
  # ============================================
  
  Pr_MHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # ORDER is important...
    cumul_surv <- exp( -Z_1 * (B2-B1))
    (B2 > B1) * cumul_surv / Nf[ B2]
  })
  
  # ===================================================
  # Paternal HSPs dynamics, survival, and probabilities
  # ===================================================
  # Population dynamics
  
  Pr_PHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # ORDER is important...
    cumul_surv <- exp( -Z_1 * (B2-B1))
    (B2 >= B1) * cumul_surv / Nm[ B2]
  })
  
  if( !all( is.finite( Pr_MHSP) & (Pr_MHSP >= 0) & (Pr_MHSP <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  if( !all( is.finite( Pr_PHSP) & (Pr_PHSP >= 0) & (Pr_PHSP <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # ===================================================
  # Calculate probability definitely a sibling
  # ===================================================
  
  Pr_MP <- autoloop(B1 = Bju_range, B2 = Bju_range, {
    # Mother father are shared
    and_mum <- alpha + log(Pr_MHSP[B1, B2] /  (1 - Pr_MHSP[B1, B2]))
    (B2 > B1) * ( Pr_PHSP[B1, B2] * (1 / (1 + exp(-and_mum))) )
  })
  
  Pr_MbarP <- autoloop(B1 = Bju_range, B2 = Bju_range, {
    # Event don't share a mother but a father
    (B2 >= B1) * (Pr_PHSP[B1, B2] - Pr_MP[B1, B2])
  })
  
  Pr_MPbar <- autoloop(B1 = Bju_range, B2 = Bju_range, {
    # Event don't share a father but a mother
    (B2 > B1) * (Pr_MHSP[B1, B2] - Pr_MP[B1, B2])
  })
  
  Pr_DSib <- autoloop(B1 = Bju_range, B2 = Bju_range, {
    # Probability of being a sibling
    (B2 >= B1) * (Pr_MP[B1, B2] + Pr_MbarP[B1, B2] * fneg_rate +  Pr_MPbar[B1, B2] * fneg_rate )
  })
  
  if( !all( is.finite( Pr_DSib) & (Pr_DSib >= 0) & (Pr_DSib <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  # -----------------
  # POP contributions
  # -----------------
  lglk <- 0
  # MOP
  lglk <- lglk + sum( dbinom(data$kin_mop, size = data$comb_mop, prob = Pr_MOP_noage, log = TRUE))
  # FOP
  lglk <- lglk + sum( dbinom(data$kin_fop, size = data$comb_fop, prob = Pr_FOP_noage, log = TRUE))
  
  # -----------------
  # HSP contributions
  # -----------------
  # SIB
  lglk <- lglk + sum(dbinom(data$kin_sib, size = data$comb_sib, prob = Pr_DSib, log = TRUE))
  
  # ------------------------------------------------
  # Refinement loop to correct for the ma/pa-ternity
  # ------------------------------------------------
  
  for (i in seq(1, dim(sib_pairs)[1]))
  {
    #print(i)
    sib_stus <- sib_pairs$FSP[i]
    bc1      <- sib_pairs$IND1_BC[i]
    bc2      <- sib_pairs$IND2_BC[i]
    pi_h2    <- sib_pairs$HAP1_FRQ[i]
    sm_hap   <- sib_pairs$SAME_HAP[i]
    # 
    # if (bc1 == bc2)
    # {
    #   next
    # }
    if (sib_stus == 1) # Full sib addition
    {
      lglk_fsib_i <- Pr_MP[bc1, bc2] / Pr_DSib[bc1, bc2]
      lglk        <- lglk + as.numeric(log(lglk_fsib_i))
     # print(paste0("Full sub contribution", log(lglk_fsib_i)))
    } else if (sib_stus == 0)
    {
      
      lglk_hsib_i <- ( (Pr_MbarP[bc1, bc2] / Pr_DSib[bc1, bc2]) * pi_h2 ) + # Paternal component
                     ( (Pr_MPbar[bc1, bc2] / Pr_DSib[bc1, bc2]) * sm_hap)  # Maternal component
      lglk        <- lglk + as.numeric(log(lglk_hsib_i))
       # if (sm_hap == 1)
       # {
       #   print(paste0("Maternal Half sub contribution", log(lglk_hsib_i)))
       # } else {
       #   print(paste0("Paternal Half sub contribution", log(lglk_hsib_i)))
       # }
      if (is.na(log(lglk_hsib_i)))
      {
        print("BAAD")
        print(pars)
      }
    }
  }
  
  # Survival prior
  if (data$surv_prior == 1)
  {
    print("HELP SUVIVAL PRIOR")
    lglk <- lglk - ((Z_1 - Mu_Z1)^2) / (2.0 * SD_Z1^2)
  }
  
  print(pars)
  print(lglk)
  if (ret_mats == 0)
  {
    return(lglk = lglk)
  } else {
    exp_mops  <- Pr_MOP_noage * data$comb_mop
    exp_fops  <- Pr_FOP_noage * data$comb_fop
    exp_fsps  <- Pr_MP * data$comb_sib
    exp_mhsps <- Pr_MPbar * data$comb_sib
    exp_fhsps <- Pr_MbarP * data$comb_sib
    exp_mats <- list(mops  = exp_mops,
                     fops  = exp_fops, 
                     fsps  = exp_fsps, 
                     mhsps = exp_mhsps, 
                     fhsps = exp_fhsps)
    
    return(list(lglk, exp_mats))
  }
}





