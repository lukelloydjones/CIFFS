 # pars <- c(500, 500, 0.3, 0.3, 0.01, 0.01)
 # years <- com.kin$years
 # y0    <- com.kin$y0
 # Bju_range <- com.kin$Bju_range
 # Yad_range <- com.kin$Yad_range
 # n_MOP_noage <- com.kin$n_MOP_noage
 # n_comp_MOP_noage <- com.kin$n_comp_MOP_noage
ciff_ckmr_loglik <- function(pars, data)
{
  N0_f  <- pars[1]
  N0_m  <- pars[2]
  Z.1   <- pars[3]
  Z.2   <- pars[4]
  R.1   <- pars[5]
  R.2   <- pars[6]
  # ============================================
  # Female dynamics, survival, and probabilities
  # ============================================
  ## Population dynamics

  Nf <- offarray( 0, first=min( years), last=max( years))
  Nf <- autoloop(yrs = years, {
    N0_f * exp( R.1 * (yrs-y0))
  })


  # offarray::autoloop() evaluates an expression over all combos of the args
  # ... and returns an array indexed by all combos. Like using nested loops.
  # Indices must match n_comp[] and n_POP[], which are set up by boring_data_prep()
  
  Pr_MOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
      cumul_surv_1 <- exp( -Z.1 * (Yad - Bju))
      cumul_surv_2 <- exp( -Z.1 * (Bju - Yad))
      if (lethal)
      {
        (Bju < Yad) * cumul_surv_1 / Nf[ Yad]
      } else {
        (Bju < Yad) * cumul_surv_1 / Nf[ Yad] + (Bju > Yad) * cumul_surv_2 / Nf[ Bju]
      }
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
  
  Nm <- offarray( 0, first=min( years), last=max( years))
  Nm <- autoloop(yrs=years, {
    N0_m * exp( R.2 * (yrs-y0))
  })

  # # offarray::autoloop() evaluates an expression over all combos of the args
  # # ... and returns an array indexed by all combos. Like using nested loops.
  # # Indices must match n_comp[] and n_POP[], which are set up by boring_data_prep()

  Pr_FOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    cumul_surv_1 <- exp( -Z.2 * (Yad - Bju))
    cumul_surv_2 <- exp( -Z.2 * (Bju - Yad))
    if (lethal)
    {
      (Bju < Yad) * cumul_surv_1 / Nm[ Yad]
    } else {
      (Bju < Yad) * cumul_surv_1 / Nm[ Yad] + (Bju > Yad) * cumul_surv_2 / Nm[ Bju]
    }
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
    (B2 > B1) * cumul_surv / Nf[ B2]
  })

  # ============================================
  # Paternal HSPs dynamics, survival, and probabilities
  # ============================================
  # Population dynamics

  Pr_PHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # ORDER is important...
    cumul_surv <- exp( -Z.2 * (B2-B1))
    (B2 > B1) * cumul_surv / Nm[ B2]
  })

  if( !all( is.finite( Pr_MHSP) & (Pr_MHSP >= 0) & (Pr_MHSP <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }
  
  if( !all( is.finite( Pr_PHSP) & (Pr_PHSP >= 0) & (Pr_PHSP <= 1))) {
    return( 1e10) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  # -----------------
  # POP contributions
  # -----------------
  lglk <- 0
  # MOP
  lglk <- lglk - sum( dbinom( n_MOP_noage, size = n_comp_MOP_noage, prob = Pr_MOP_noage, log = TRUE))
  # FOP
  lglk <- lglk - sum( dbinom( n_FOP_noage, size = n_comp_FOP_noage, prob = Pr_FOP_noage, log = TRUE))
 
  # -----------------
  # HSP contributions
  # -----------------
  # MHSP
  lglk <- lglk - sum( dbinom( n_mhsps, size = n_comp_HSP, prob = Pr_MHSP, log = TRUE))
  # PHSP
  lglk <- lglk - sum( dbinom( n_phsps, size = n_comp_HSP, prob = Pr_PHSP, log = TRUE))
  
  N <- Nf + Nm
  list2env( mget( cq( Nf, Nm, N, Pr_MOP_noage, Pr_FOP_noage, Pr_MHSP, Pr_PHSP)), com.kin)
  return(lglk)
}