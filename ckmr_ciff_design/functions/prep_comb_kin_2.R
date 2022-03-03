# -----------------------------------------------------
# Function to compute combinations and kin numbers 
# for each year and kin type
# More Markovian. THIS IS VERY IMPORTANT TO GET CORRECT
# -----------------------------------------------------
prep_comb_kin_2 <- function(samp.no.py, years.samp, prop_ju_py, prop_sa_py, prop_ad_py, f_prop, lethal)
{
  library(dplyr)
  # PARS
  #  kin <- the set of indices of the sample that are of a specific kin type
  #  pop.samp - population sample for rep i
  #  lethal - prep data as if they are lethal mum's sampling year > birth of offspirng
  #           if 0 the non-lethal birth > year as she can be sampled before
  # RETURNS
  #  comb.kin - four dimensional array woth combinations, no. kin 
  #             and types
  # -------------------------------------------------------------
  # Who are juveniles, who are juveniles or adults, who are just adults
  # ------------------------
  end_yr = 2018
  m_prop = 1-f_prop

  # -----------------------------------------------------------
  # Only compare known adults with known juveniles or subadults
  # -----------------------------------------------------------
  
  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Yad_range <- (end_yr - years.samp + 1)  %upto% end_yr
  Bju_range <- (min(Yad_range)-2) %upto% max((Yad_range)-1)
  
  # Package up stuff, to be used as environment for lglk function
  y0    <- min( Bju_range)
  years <- min( Bju_range) %upto% max( Yad_range)
  
  # Construct 'sampling' from samples per year + proportions
  
  # mf_... is samp size. FEMALES
  mf_ju_B  <- offarray(c(0, round(samp.no.py * prop_ju_py * f_prop)   ), first=min(Bju_range), last=max(Bju_range))
  mf_sa_B  <- offarray(c(   round(samp.no.py * prop_sa_py * f_prop), 0), first=min(Bju_range), last=max(Bju_range))
  mf_ad_Y  <- offarray(c(   round(samp.no.py * prop_ad_py * f_prop)   ), first=min(Yad_range), last=max(Yad_range))
  
  # mf_... is samp size. MALES
  mm_ju_B  <- offarray(c(0, round(samp.no.py * prop_ju_py * m_prop)   ), first=min(Bju_range), last=max(Bju_range))
  mm_sa_B  <- offarray(c(   round(samp.no.py * prop_sa_py * m_prop), 0), first=min(Bju_range), last=max(Bju_range))
  mm_ad_Y  <- offarray(c(   round(samp.no.py * prop_ad_py * m_prop)   ), first=min(Yad_range), last=max(Yad_range))
  
  # m_... is samp size. GLOBAL
  m_ju_B    <- mf_ju_B + mm_ju_B
  m_sa_B    <- mf_sa_B + mm_sa_B
  m_ju_sa_B <- m_ju_B  + m_sa_B
  m_ad_Y    <- mf_ad_Y + mm_ad_Y
  
  # ------
  # Global - not used. Just for checking
  # ------
  
  # n_comp_POP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
  #   # No point if ju born after adult was (lethally) sampled
  #   m_ju_sa_B[ Bju] * m_ad_Y[ Yad] * (Bju <  Yad)
  # })
  
  # --------------------------------------
  # MOPS
  # --------------------------------------
  n_comp_MOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    # No point if ju born after adult was (lethally) sampled
    if (lethal)
    {
      m_ju_sa_B[ Bju] * mf_ad_Y[ Yad] * (Bju < Yad)
    } else {
      m_ju_sa_B[ Bju] * (mf_ad_Y[ Yad] + (mf_sa_B[Yad-2] * (Bju - Yad >= 1))  + (mf_ju_B[Yad-1] * (Bju - Yad >= 2))) * (Bju != Yad)
    }
  })
  
  # --------------------------------------
  # FOPS
  # --------------------------------------
  n_comp_FOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    # No point if ju born after adult was (lethally) sampled
    if (lethal)
    {
      m_ju_sa_B[ Bju] * mm_ad_Y[ Yad] * (Bju <  Yad)
    } else {
      m_ju_sa_B[ Bju] * (mm_ad_Y[ Yad] + (mm_sa_B[Yad-2] * (Bju - Yad >= 1))  + (mm_ju_B[Yad-1] * (Bju - Yad >= 2))) * (Bju != Yad)
    }
  })
  
  # n_comp_MOP_noage + n_comp_FOP_noage

  # --------------------------------------
  # Half-siblings
  # --------------------------------------
  
  # And for HSPs... guaranteed in birth-order by row, thx2 prepare_from_sim()
  n_comp_HSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # NB *exclude* double-count and same-cohort
    m_ju_sa_B[ B1] * m_ju_sa_B[ B2] * (B2>B1)
  })


  envo <- list2env( mget( cq(y0, years)), parent=.GlobalEnv)
  list2env( mget( cq(n_comp_MOP_noage,
                     n_comp_FOP_noage,
                     n_comp_HSP, 
                     Bju_range, Yad_range, lethal)), envo)
  return(envo)
}
