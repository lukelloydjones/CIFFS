#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(comb_mop);                         // Data vector transmitted from R
  DATA_MATRIX(comb_fop);                         // Data vector transmitted from R
  DATA_MATRIX(comb_sib);                         // combination matrix for siblings at broad scale

    
  DATA_MATRIX(kin_mop);                          // Data vector transmitted from R
  DATA_MATRIX(kin_fop);                          // Data vector transmitted from R
  DATA_MATRIX(kin_sib);                          // Sibling kin matrix at the broad scale
  
  // Data vectors on siblings for mtDNA clean-up loop.

  DATA_IVECTOR(dsm_bc1);                         // Birth cohort IND1 as an integer from bymin
  DATA_IVECTOR(dsm_bc2);                         // Birth cohort IND2 as an integer from bymin
  DATA_VECTOR(dsm_pi_h2);                        // Frequency of mtDNA haplotype for IND2                   
  DATA_IVECTOR(dsm_sm_hp);                       // Bool (0/1) for whether pair share same haplotype. 1 if yes                          
  DATA_IVECTOR(df_fsp_yn);                       // Bool (0/1) for whether pair are FSP or HSP. 1 if FSP.                                                                         

  DATA_INTEGER(symin);                           // Data vector transmitted from R
  DATA_INTEGER(symax);                           // Data vector transmitted from R
  DATA_INTEGER(bymin);                           // Data vector transmitted from R
  DATA_INTEGER(bymax);                           // Data vector transmitted from R
  DATA_INTEGER(lethal);                          // Data vector transmitted from R
  DATA_INTEGER(surv_prior);                      // Data vector transmitted from R
  DATA_INTEGER(roi_prior);                       // Data vector transmitted from R
    
  DATA_SCALAR(fng_rate);                         // False negative rate from the threshold of PLODogram

  DATA_SCALAR(Mu_Z1);                            // Normal prior mean for survival if required
  DATA_SCALAR(SD_Z1);                            // Normal prior mean for survival if required

  DATA_SCALAR(Mu_R1);                            // Normal prior mean for ROI if required
  DATA_SCALAR(SD_R1);                            // Normal prior variance for ROI if required
    
  PARAMETER(N0_f);                               // Females abundance year naught
  PARAMETER(IP);                                 // Male abundance supression factor
  PARAMETER(Z1);                                 // Female/or and male survival
  //PARAMETER(Z2);                                 // Female and male survival
  PARAMETER(R1);                                 // Female/or and male ROC
  //PARAMETER(R2);                                 // Male ROC
  PARAMETER(Alpha);                              // Assortativity parameter


  Type ciff_llk;                                 // TMB weird declare the "objective function" (neg. log. likelihood)

  // Initial declarations
  int yrs_tot, byrs_tot, syrs_tot, year, byr_1, byr_2, syrs_diff, yad, bju, rel_yad, sib_mt_it;
  int df_sib_tot = dsm_pi_h2.size();
  yrs_tot   = symax - bymin + 1;
  syrs_tot  = symax - symin + 1;
  byrs_tot  = bymax - bymin + 1;
  syrs_diff = symin - bymin;

  std::cout << "Years total " << yrs_tot << std::endl;
  std::cout << "syrs_tot "    << syrs_tot << std::endl;
  std::cout << "byrs_tot "     << byrs_tot << std::endl;
  std::cout << "Number of deffo sibs " << df_sib_tot << std::endl;
    
  vector<Type> N_f(yrs_tot);
  vector<Type> N_m(yrs_tot);
    
  matrix<Type> mop_prob(byrs_tot,  syrs_tot);
  matrix<Type> exp_mops(byrs_tot,  syrs_tot);
  matrix<Type> fop_prob(byrs_tot,  syrs_tot);
  matrix<Type> mhsp_prob(byrs_tot, byrs_tot);
  matrix<Type> m_p_prob(byrs_tot,  byrs_tot);
  matrix<Type> mbar_p_prob(byrs_tot,  byrs_tot);
  matrix<Type> m_pbar_prob(byrs_tot,  byrs_tot);
  matrix<Type> df_sib_prob(byrs_tot,  byrs_tot);

  matrix<Type> exp_mhsps(byrs_tot, byrs_tot);
  matrix<Type> fhsp_prob(byrs_tot, byrs_tot);

  // Year total span from the first birth year to the last sample year e.g., samples years 2015-2018, birth years 2010,2017
  // The cumulative survival is computed for years 2010-2018. Remember 0 index
  // The index for the sampling years is yrs_tot - syrs_tot  to syrs_tot. Check yrs_tot = (2018 - 2010) = 8
  // The index for the birth years is yrs_tot - syrs_tot  to syrs_tot. Check yrs_tot = (2018 - 2010) = 8

  // ===============
  // Female dynamics
  // ===============
  // R1 = exp(R1); Change of scale for optimisation?
  N_f[0] = exp( N0_f );
  for (year = 1; year < yrs_tot; year++)
  {
    N_f[year] = N_f[0] * exp(R1 * year);
  }
    
  // =============
  // Male dynamics 
  // =============
   N_m[0] = exp(N0_f) * (exp(IP) / (1 + exp(IP)));
   for (year = 1; year < yrs_tot; year++)
   {
     N_m[year] = N_m[0] * exp(R1 * year);
     //N_m[year] = N_m[0] * exp(R2 * year);
   }

  // ======================
  // CKMR POP probabilities
  // ======================

   for (yad = 0; yad < syrs_tot; yad++)
   {
     rel_yad = yad + syrs_diff;
     for (bju = 0; bju < byrs_tot; bju++)
     {
       // Female survival
       Type cumul_surv_f1 = exp(-Z1 * (rel_yad - bju));
       Type cumul_surv_f2 = exp(-Z1 * (bju - rel_yad));
       // Male survival
       Type cumul_surv_m1 = exp(-Z1 * (rel_yad - bju));
       Type cumul_surv_m2 = exp(-Z1 * (bju - rel_yad));
      
       //std::cout << "Proposed index " << rel_yad << std::endl;
       //std::cout << "Indices " << yad  << " " << bju << std::endl;
       if (bju < rel_yad)
       {
         //std::cout << "POP prob 1" << std::endl;
         mop_prob(bju, yad) = cumul_surv_f1 / N_f[rel_yad];
         fop_prob(bju, yad) = cumul_surv_m1 / N_m[rel_yad];
       } else if (!lethal && (bju > rel_yad))
       {
         //std::cout << "POP prob 2" << std::endl;
         mop_prob(bju, yad) = cumul_surv_f2 / N_f[bju];
         fop_prob(bju, yad) = cumul_surv_m2 / N_m[bju];
       }
     }
  }

   std::cout << "MOP probabilities" << std::endl;
   for (year = 0; year < byrs_tot; year++)
   {
     std::cout << mop_prob.row(year) << std::endl;
   }

   std::cout << "FOP probabilities" << std::endl;
   for (year = 0; year < byrs_tot; year++)
   {
     std::cout << fop_prob.row(year) << std::endl;
   }
  
  // ==================
  // Log-likelihood POP
  // ==================
  //  std::cout << "Starting POP llk" << std::endl;
  //  std::cout << "Sampling years total"  << syrs_tot << std::endl;
  //  std::cout << "Birth years total" << byrs_tot<< std::endl;
  ciff_llk = 0.0;
  for (yad = 0; yad < syrs_tot; yad++)
  {
	  rel_yad = yad + syrs_diff;
    for (bju = 0; bju < byrs_tot; bju++)
    {
	    if (bju < rel_yad || (!lethal && bju > rel_yad))
	    {
	 	   ciff_llk -= dbinom(kin_mop(bju, yad),
	 		   			          	comb_mop(bju, yad),
	 		 			          	  mop_prob(bju, yad), true);
      // exp_mops(bju, yad) = comb_mop(bju, yad) * mop_prob(bju, yad);
	 	  //std::cout << "MOP llk " << ciff_llk << std::endl;
	 	  //std::cout << "year birth year " << yad << " " << bju <<
	 		//		 " kin " << kin_mop(bju, yad) <<
	 		//		 " combs " << comb_mop(bju, yad) <<
	 		//		 " prob " << mop_prob(bju, yad) << std::endl;
	 	    ciff_llk -= dbinom(kin_fop(bju, yad),
	 			  			           comb_fop(bju, yad),
	 			  			           fop_prob(bju, yad), true);
	    //      std::cout << "FOP llk " << ciff_llk << std::endl;
	    //      std::cout << "year birth year " << yad << " " << bju << std::endl;
	    }
    }
  }
  std::cout << "POP llk " << ciff_llk << std::endl;

  // ===============================
  // CKMR HSP maternal probabilities
  // ===============================
  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      if ((byr_1 < byr_2))
      {
        // Female survival
        Type cumul_surv_f = exp(-Z1 * (byr_2 - byr_1)); 
        mhsp_prob(byr_1, byr_2) = cumul_surv_f / N_f[byr_2];
      }
    }
  }

  // ===============================
  // CKMR HSP paternal probabilities
  // ===============================
  
  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      if ((byr_1 <= byr_2))
      {
        // Male survival
        Type cumul_surv_m = exp(-Z1 * (byr_2 - byr_1));
        fhsp_prob(byr_1, byr_2) = cumul_surv_m / N_m[byr_2];
      }
    }
  }

  // ========================================
  // Compute the deffo sibling probabilities
  // ========================================

  // Maternal AND Paternal

  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      if ((byr_1 <= byr_2))
      {
        // maternal and paternal
        Type and_mum = Alpha + log(mhsp_prob(byr_1, byr_2) / (1.0 - mhsp_prob(byr_1, byr_2))); // Assortativity component
        m_p_prob(byr_1, byr_2) = fhsp_prob(byr_1, byr_2) * (1 / (1 + exp(-and_mum)));
      }
    }
  }

  // Not Maternal AND Paternal

  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      if ((byr_1 <= byr_2))
      {
        mbar_p_prob(byr_1, byr_2) = fhsp_prob(byr_1, byr_2) - m_p_prob(byr_1, byr_2);
      }
    }
  }

  // Maternal AND Not Paternal

  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      if ((byr_1 < byr_2))
      {
        m_pbar_prob(byr_1, byr_2) = mhsp_prob(byr_1, byr_2) - m_p_prob(byr_1, byr_2);
      }
    }
  }

  // Definitely a sibling probability

  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      if ((byr_1 <= byr_2))
      {
        df_sib_prob(byr_1, byr_2) =  ( m_p_prob(byr_1, byr_2) + 
                                       mbar_p_prob(byr_1, byr_2) * fng_rate +  
                                       m_pbar_prob(byr_1, byr_2) * fng_rate ) ;
      }
    }
  }

  std::cout << "Def sib probabilities" << std::endl;
  for (year = 0; year < byrs_tot; year++)
  {
     std::cout << df_sib_prob.row(year) << std::endl;
  }


  // ================================
  // Log-likelihood HSP - Broad scale
  // ================================
  std::cout << "Starting HSP llk" << std::endl;
  for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
    {
      // if (byr_1 <= byr_2)
      if (byr_1 < byr_2)
      {
        // std::cout << "HSP llk " << ciff_llk << std::endl;
        // std::cout << "by1 by2 " << byr_1 << " " << byr_2 << std::endl;
        ciff_llk -= dbinom( kin_sib(byr_1, byr_2),
                            comb_sib(byr_1, byr_2),
                            df_sib_prob(byr_1, byr_2), true);
      }
    }
  }

  std::cout << "Total llk " << ciff_llk << std::endl;
  // =============================================
  // Have another loop here for the mtDNA clean up
  // =============================================
  // Declare some variables. Should do it at top.
  int  bc1, bc2, sm_hap, sib_stus;
  Type pi_h2, lglk_fsib_i, lglk_hsib_i;                                                                      

  for (sib_mt_it = 0; sib_mt_it < df_sib_tot; sib_mt_it++)
  {
    
    bc1      = dsm_bc1(sib_mt_it);
    bc2      = dsm_bc2(sib_mt_it);
    pi_h2    = dsm_pi_h2(sib_mt_it);
    sm_hap   = dsm_sm_hp(sib_mt_it);
    sib_stus = df_fsp_yn(sib_mt_it);
    
    if (sib_stus == 1) // Full sib contribution
    {

      lglk_fsib_i = m_p_prob(bc1, bc2) / df_sib_prob(bc1, bc2);
      ciff_llk    = ciff_llk - log(lglk_fsib_i);

      std::cout << "Full sib contribution " << log(lglk_fsib_i) << std::endl;

    } else if (sib_stus == 0) { // Half sib contribution

      lglk_hsib_i = ( (mbar_p_prob(bc1, bc2) / df_sib_prob(bc1, bc2) * pi_h2 ) +  // Paternal component
                      (m_pbar_prob(bc1, bc2) / df_sib_prob(bc1, bc2) * sm_hap) ); // Maternal component
      ciff_llk    = ciff_llk - log(lglk_hsib_i);

      if (sm_hap == 1)
      {
        std::cout << "Maternal Half sib contribution " << log(lglk_hsib_i) << std::endl;
      } else {
        std::cout << "Paternal Half sib contribution " << log(lglk_hsib_i) << std::endl;
      }

    }
  }
  // ======================
  // Survival prior
  // ======================
  if (surv_prior)
  {
    std::cout << "HEEELLP SURVIVAL PRIOR" << std::endl;
    ciff_llk += pow(Z1 - Mu_Z1, 2.0) / (2.0 * pow(SD_Z1, 2.0));
    // ciff_llk += pow(Z2 - Mu_Z1, 2.0) / (2.0 * pow(SD_Z1, 2.0));
  }


  if (roi_prior)
  {
    std::cout << "HEEELLP ROI PRIOR" << std::endl;
    ciff_llk += pow(R1 - Mu_R1, 2.0) / (2.0 * pow(SD_R1, 2.0));
    // ciff_llk += pow(Z2 - Mu_Z1, 2.0) / (2.0 * pow(SD_Z1, 2.0));
  }

  REPORT(mop_prob);
  REPORT(mhsp_prob);
  REPORT(exp_mops);
  REPORT(exp_mhsps);
  ADREPORT(N_f);
  ADREPORT(log(N_f));
  ADREPORT(N_m);
  ADREPORT(1 / (1 + exp(-IP)))
  //ADREPORT(log(N_m));
  //ADREPORT(log(N_f + N_m));
  return ciff_llk;
}
