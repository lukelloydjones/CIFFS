# -----------------------------------------------------
# Function to compute combinations and kin numbers 
# for each year and kin type
# More Markovian. THIS IS VERY IMPORTANT TO GET CORRECT
# -----------------------------------------------------
prep_comb_kin <- function(kin, pop.samp.i, lethal)
{
  library(dplyr)
  # PARS
  #  kin <- the set of indices of the sample that are of a specific kin type
  #  pop.samp.i - population sample for rep i
  #  lethal - prep data as if they are lethal mum's sampling year > birth of offspirng
  #           if 0 the non-lethal birth > year as she can be sampled before
  # RETURNS
  #  comb.kin - four dimensional array woth combinations, no. kin 
  #             and types
  # -------------------------------------------------------------
  # Who are juveniles, who are juveniles or adults, who are just adults
  # Lethal sampling scenario
  # ------------------------
  Y <- pop.samp.i$year_samp
  #A <- pop.samp.i$ages
  B <- pop.samp.i$year #-1 #- 2 
  
  # order the pairs
  mops  <- kin$pops
  #fops  <- kin$fops
  mhsps <- kin$hsps
  #phsps <- kin$phsps
  
  dim(mops)
  #dim(fops)
  dim(mhsps)
  #dim(phsps)
  
  # -----------------------------------------------------------
  # Only compare known adults with known juveniles or subadults
  # -----------------------------------------------------------
  juvid     <- which(pop.samp.i$life.sg.status == "J")
  juvid.f   <- which(pop.samp.i$life.sg.status == "J" & pop.samp.i$mf == 0)
  juvid.m   <- which(pop.samp.i$life.sg.status == "J" & pop.samp.i$mf == 1)
  subadid   <- which(pop.samp.i$life.sg.status == "SA")
  subadid.f <- which(pop.samp.i$life.sg.status == "SA" & pop.samp.i$mf == 0)
  subadid.m <- which(pop.samp.i$life.sg.status == "SA" & pop.samp.i$mf == 1)
  adid      <- which(pop.samp.i$life.sg.status == "A")
  adid.f    <- which(pop.samp.i$life.sg.status == "A" & pop.samp.i$mf == 0)
  adid.m    <- which(pop.samp.i$life.sg.status == "A" & pop.samp.i$mf == 1)
  juv.subad.id <- which(pop.samp.i$life.sg.status == "J" | pop.samp.i$life.sg.status == "SA")

  # Put in birth order
  mops[(B[ mops[, 1]] > B[ mops[, 2]]), ] <- cbind(mops[(B[ mops[, 1]] > B[ mops[, 2]]), 2], mops[(B[ mops[, 1]] > B[ mops[, 2]]), 1])
  #fops[(B[ fops[, 1]] > B[ fops[, 2]]), ] <- cbind(fops[(B[ fops[, 1]] > B[ fops[, 2]]), 2], fops[(B[ fops[, 1]] > B[ fops[, 2]]), 1])
  
  # Remove adult-adult comparisons. 06/12/2019 we moved this check to the find kin function
  if (length(which(paste0(pop.samp.i$life.sg.status[mops[, 1]], pop.samp.i$life.sg.status[mops[, 2]]) == "AA"))> 0)
  {
    mops <- mops[-which(paste0(pop.samp.i$life.sg.status[mops[, 1]], pop.samp.i$life.sg.status[mops[, 2]]) == "AA"), ]
  }
    #fops <- fops[-which(paste0(pop.samp.i$life.sg.status[fops[, 1]], pop.samp.i$life.sg.status[fops[, 2]]) == "AA"), ]
  
  dim(mops)
  # dim(fops)
  
  #mops.pops <- cbind(pop.samp.i[mops[, 1],], pop.samp.i[mops[, 2], ])
  #head(mops.pops[(mops.pops[, 8] < mops.pops[, 15]), ], 50)
  
  if (lethal)
  {
    mops <- mops[ Y[ mops[,1]] > B[ mops[, 2]],]
    #fops <- fops[ Y[ fops[,1]] > B[ fops[, 2]],]
  } else if (!lethal)
  {
    if (length(which(pop.samp.i$life.sg.status[mops[, 2]] == "A")) > 0)
    {
      mops <- mops[-which(pop.samp.i$life.sg.status[mops[, 2]] == "A"), ]
    } 
    #else if (length(which(pop.samp.i$life.sg.status[fops[, 2]] == "A")) > 0)
    #{
      #fops <- fops[-which(pop.samp.i$life.sg.status[fops[, 2]] == "A"), ]
    #}
  }

  # dim(mops)
  # dim(fops)
  # 06/12/2019. Done again in the find_kin function. We hope
  # if (!lethal)
  # {
  #   mops <- mops[-which(duplicated(mops)), ]
  #   fops <- fops[-which(duplicated(fops)), ]
  # }
  
  # ----------
  # HSPs
  # ----------
  # Put in birth order
  if (sum((B[ mhsps[, 1]] > B[ mhsps[, 2]])) > 0)
  {
    mhsps[(B[ mhsps[, 1]] > B[ mhsps[, 2]]), ] <- cbind(mhsps[(B[ mhsps[, 1]] > B[ mhsps[, 2]]), 2], mhsps[(B[ mhsps[, 1]] > B[ mhsps[, 2]]), 1])
  }
  # if (sum((B[ phsps[, 1]] > B[ phsps[, 2]])) > 0)
  # {
  #   phsps[(B[ phsps[, 1]] > B[ phsps[, 2]]), ] <- cbind(phsps[(B[ phsps[, 1]] > B[ phsps[, 2]]), 2], phsps[(B[ phsps[, 1]] > B[ phsps[, 2]]), 1])
  # }
  
  mhsps <- mhsps[which(paste0(pop.samp.i$life.sg.status[mhsps[, 1]], pop.samp.i$life.sg.status[mhsps[, 2]]) == "JJ"  | 
                       paste0(pop.samp.i$life.sg.status[mhsps[, 1]], pop.samp.i$life.sg.status[mhsps[, 2]]) == "JSA" |
                       paste0(pop.samp.i$life.sg.status[mhsps[, 1]], pop.samp.i$life.sg.status[mhsps[, 2]]) == "SAJ" |
                       paste0(pop.samp.i$life.sg.status[mhsps[, 1]], pop.samp.i$life.sg.status[mhsps[, 2]]) == "SASA"), ]
  # phsps <- phsps[which(paste0(pop.samp.i$life.sg.status[phsps[, 1]], pop.samp.i$life.sg.status[phsps[, 2]]) == "JJ"  | 
  #                      paste0(pop.samp.i$life.sg.status[phsps[, 1]], pop.samp.i$life.sg.status[phsps[, 2]]) == "JSA" |
  #                      paste0(pop.samp.i$life.sg.status[phsps[, 1]], pop.samp.i$life.sg.status[phsps[, 2]]) == "SAJ" |
  #                      paste0(pop.samp.i$life.sg.status[phsps[, 1]], pop.samp.i$life.sg.status[phsps[, 2]]) == "SASA"), ]


  print(dim(mops))
  #dim(fops)
  print(dim(mhsps))
  #dim(phsps)
  
  #mhsps <- mhsps[-which(duplicated(mhsps)), ]
  #phsps <- phsps[-which(duplicated(phsps)), ]
  
  #dim(mops)
  #dim(fops)
  dim(mhsps)
  #dim(phsps)
  
  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Yad_range <- min( Y[ adid])  %upto% max( Y[ adid])
  #Aad_range <- min( A[ adid])  %upto% max( A[ adid])
  Bju_range <- (min(Yad_range)-2) %upto% max( B[ juv.subad.id]) # In case no subads are sampled in year 1
  #Bju_range <- min(B[ juv.subad.id]) %upto% max( B[ juv.subad.id]) # In case no subads are sampled in year 1 # POTENTIAL CHANGE ON 4/02/2020
  
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
  
  n_comp_POP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    # No point if ju born after adult was (lethally) sampled
    m_ju_sa_B[ Bju] * m_ad_Y[ Yad] * (Bju <  Yad)
  })
  
  # # --------------------------------------
  # # MOPS
  # # --------------------------------------
  n_comp_MOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    # No point if ju born after adult was (lethally) sampled
    if (lethal)
    {
      m_ju_sa_B[ Bju] * mf_ad_Y[ Yad] * (Bju < Yad)
    } else {
      m_ju_sa_B[ Bju] * (mf_ad_Y[ Yad] + (mf_sa_B[Yad-2] * (Bju - Yad >= 1))  + (mf_ju_B[Yad-1] * (Bju - Yad >= 2))) * (Bju != Yad)
    }
  })
  # # 
  n_MOP_noage_all <- offarray( table( Bju=B[ mops[,2]], Yad=Y[ mops[,1]]),
                              template=n_comp_MOP_noage)
  # Lethal-non lethal clean up
  n_MOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    # No point if ju born after adult was (lethally) sampled
    if (lethal)
    {
      n_MOP_noage_all[ Bju, Yad]  * (Bju < Yad)
    } else {
      n_MOP_noage_all[ Bju, Yad] * (Bju != Yad)
    }
  })
  # # --------------------------------------
  # # FOPS
  # # --------------------------------------
  n_comp_FOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
    # No point if ju born after adult was (lethally) sampled
    if (lethal)
    {
      m_ju_sa_B[ Bju] * mm_ad_Y[ Yad] * (Bju <  Yad)
    } else {
      m_ju_sa_B[ Bju] * (mm_ad_Y[ Yad] + (mm_sa_B[Yad-2] * (Bju - Yad >= 1))  + (mm_ju_B[Yad-1] * (Bju - Yad >= 2))) * (Bju != Yad)
    }
  })
  # n_comp_MOP_noage <- n_comp_MOP_noage + n_comp_FOP_noage # This is a fudge thing when just look at not seprating by sex # Date 04/02/2020
  # 
  # # Now total up number of POPs seen, in the same way
  # n_FOP_noage_all <- offarray( table( Bju=B[ fops[,2]], Yad=Y[ fops[,1]]),
  #                          template=n_comp_FOP_noage)
  # 
  # # Lethal-non lethal clean up
  # n_FOP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
  #   # No point if ju born after adult was (lethally) sampled
  #   if (lethal)
  #   {
  #     n_FOP_noage_all[ Bju, Yad]  * (Bju < Yad)
  #   } else {
  #     n_FOP_noage_all[ Bju, Yad] * (Bju != Yad) 
  #   }
  # })
  # --------------------------------------
  # Half-siblings
  # --------------------------------------
  
  # And for HSPs... guaranteed in birth-order by row, thx2 prepare_from_sim()
  n_comp_HSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    # NB *exclude* double-count and same-cohort
    m_ju_sa_B[ B1] * m_ju_sa_B[ B2] * (B2>=B1)
  })
  
  n_mhsps <- offarray( table( B1=B[ mhsps[,1]], B2=B[ mhsps[,2]]),
                     template=n_comp_HSP)

  
  #n_phsps <- offarray( table( B1=B[ phsps[,1]], B2=B[ phsps[,2]]),
  #                   template=n_comp_HSP)

  envo <- list2env( mget( cq(B, Y, y0, years, juv.subad.id, adid)), parent=.GlobalEnv)
  list2env( mget( cq(n_comp_MOP_noage, n_comp_FOP_noage, n_MOP_noage, n_comp_HSP, n_mhsps, 
                     Bju_range, Yad_range, lethal)), envo)
  return(envo)
}
