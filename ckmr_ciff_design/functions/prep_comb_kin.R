# -----------------------------------------------------
# Function to compute combinations and kin numbers 
# for each year and kin type
# More Markovian. THIS IS VERY IMPORTANT TO GET CORRECT
# -----------------------------------------------------
prep_comb_kin <- function(kin, pop.samp, lethal)
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
  # Lethal sampling scenario
  # ------------------------
  Y <- pop.samp$year_samp
  A <- pop.samp$ages
  B <- pop.samp$year - 1 

  
  # -----------------------------------------------------------
  # Only compare known adults with known juveniles or subadults
  # -----------------------------------------------------------
  juvid     <- which(pop.samp$life.sg.status == "J")
  juvid.f   <- which(pop.samp$life.sg.status == "J" & pop.samp$mf == 0)
  juvid.m   <- which(pop.samp$life.sg.status == "J" & pop.samp$mf == 1)
  subadid   <- which(pop.samp$life.sg.status == "SA")
  subadid.f <- which(pop.samp$life.sg.status == "SA" & pop.samp$mf == 0)
  subadid.m <- which(pop.samp$life.sg.status == "SA" & pop.samp$mf == 1)
  adid      <- which(pop.samp$life.sg.status == "A")
  adid.f    <- which(pop.samp$life.sg.status == "A" & pop.samp$mf == 0)
  adid.m    <- which(pop.samp$life.sg.status == "A" & pop.samp$mf == 1)
  juv.subad.id <- which(pop.samp$life.sg.status == "J" | pop.samp$life.sg.status == "SA")
  
  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Aad_range <- min( A[ adid])  %upto% max( A[ adid])
  Yad_range <- min( Y[ adid])  %upto% max( Y[ adid])
  #Bju_range <- min( B[ juv.subad.id]) %upto% max( B[ juv.subad.id])
  Bju_range <- (min(Yad_range)-2) %upto% max( B[ juv.subad.id]) # In case no subads are sampled in year 1
  
  # Package up stuff, to be used as environment for lglk function
  # y0    <- min( B[ juv.subad.id]) # SHOULDN'T really be data-driven
  # years <- min( B[ juv.subad.id]) %upto% max( Y[ adid])
  y0    <- min( Bju_range) # SHOULDN'T really be data-driven
  years <- min( Bju_range) %upto% max( Yad_range)
  # #envo <- list2env( mget( cq(B, Y, A, y0, years, juv.subad.id, adid)), parent=.GlobalEnv)
  
  # Add psuedo counts to ensure dimensions, then remove after offarray construction
  psuedo_Bju = offarray(table(Bju_range))
  psuedo_Yad = offarray(table(Yad_range))
  
  # Create with psudeo counts
  
  # m_... is samp size. GLOBAL
  m_ju_sa_B <- offarray( table( c(B[ juv.subad.id], Bju_range) ) )
  m_ju_B    <- offarray( table( c(B[ juvid], Bju_range) ) )
  m_sa_B    <- offarray( table( c(B[ subadid], Bju_range) ) )
  m_ad_Y    <- offarray( table( c(Y[ adid], Yad_range) ) ) # for noage
  # m_... is samp size. FEMALES
  mf_ju_B  <- offarray( table( c(B[ juvid.f], Bju_range) ) )
  mf_sa_B  <- offarray( table( c(B[ subadid.f], Bju_range) ) )
  mf_ad_Y  <- offarray( table( c(Y[ adid.f], Yad_range) ) ) # for noage
  # m_... is samp size. MALES
  mm_ju_B  <- offarray( table( c(B[ juvid.m], Bju_range) ) )
  mm_sa_B  <- offarray( table( c(B[ subadid.m], Bju_range) ) )
  mm_ad_Y  <- offarray( table( c(Y[ adid.m], Yad_range) ) ) # for noage
  
  # Remove psuedo counts before use
  
  # m_... is samp size. GLOBAL
  m_ju_sa_B <- m_ju_sa_B - psuedo_Bju
  m_ju_B    <- m_ju_B    - psuedo_Bju
  m_sa_B    <- m_sa_B    - psuedo_Bju
  m_ad_Y    <- m_ad_Y    - psuedo_Yad # for noage
  # m_... is samp size. FEMALES
  mf_ju_B  <- mf_ju_B - psuedo_Bju
  mf_sa_B  <- mf_sa_B - psuedo_Bju
  mf_ad_Y  <- mf_ad_Y - psuedo_Yad # for noage
  # m_... is samp size. MALES
  mm_ju_B  <- mm_ju_B - psuedo_Bju
  mm_sa_B  <- mm_sa_B - psuedo_Bju
  mm_ad_Y  <- mm_ad_Y - psuedo_Yad # for noage
  
  
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


  envo <- list2env( mget( cq(B, Y, A, y0, years, juv.subad.id, adid)), parent=.GlobalEnv)
  list2env( mget( cq(n_comp_MOP_noage,
                     n_comp_FOP_noage,
                     n_comp_HSP, 
                     Bju_range, Yad_range, Aad_range, lethal)), envo)
  return(envo)
}
