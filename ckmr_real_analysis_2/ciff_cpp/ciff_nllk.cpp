#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(comb_mop);                          // Data vector transmitted from R
  //DATA_MATRIX(comb_fop);                          // Data vector transmitted from R
  DATA_MATRIX(comb_mhsp);                         // Data vector transmitted from R
  //DATA_MATRIX(comb_fhsp);                         // Data vector transmitted from R
    
  DATA_MATRIX(kin_mop);                           // Data vector transmitted from R
  //DATA_MATRIX(kin_fop);                           // Data vector transmitted from R
  DATA_MATRIX(kin_mhsp);                          // Data vector transmitted from R
  //DATA_MATRIX(kin_fhsp);                          // Data vector transmitted from R
    
  DATA_INTEGER(symin);                             // Data vector transmitted from R
  DATA_INTEGER(symax);                             // Data vector transmitted from R
  DATA_INTEGER(bymin);                            // Data vector transmitted from R
  DATA_INTEGER(bymax);                            // Data vector transmitted from R
  DATA_INTEGER(lethal);                           // Data vector transmitted from R
    
  DATA_SCALAR(fudge_factor);                           //
    
  PARAMETER(N0_f);                                // Parameter value transmitted from R - Initial number of females
  //PARAMETER(N0_m);                                // - Initial number of males
  PARAMETER(Z1);                                  // - Female mortality
  //PARAMETER(Z2);                                  // - Male mortality
  PARAMETER(R1);                                  // - Female recruitment
 // PARAMETER(R2);                                  // - Male recruitment

  Type ciff_llk;                                  // TMB weird declare the "objective function" (neg. log. likelihood)

  // Initial declarations
  int yrs_tot, byrs_tot, syrs_tot, year, birth_year_1, birth_year_2, start_yrs_diff;
  yrs_tot = symax - bymin + 1;
  syrs_tot = symax -  symin + 1;
  byrs_tot = bymax - bymin + 1;
  start_yrs_diff = symin - bymin;

  std::cout << "Years total " << yrs_tot << std::endl;
  std::cout << "syrs_tot " << syrs_tot << std::endl;
  std::cout << "byrs_tot" << byrs_tot << std::endl;
    
  vector<Type> N_f(yrs_tot);
  //vector<Type> N_m(yrs_tot);
    
  matrix<Type> mop_prob(byrs_tot, syrs_tot);
  matrix<Type> exp_mops(byrs_tot, syrs_tot);
  //matrix<Type> fop_prob(byrs_tot, syrs_tot);
  matrix<Type> mhsp_prob(byrs_tot, byrs_tot);
  matrix<Type> exp_mhsps(byrs_tot, byrs_tot);
  //matrix<Type> fhsp_prob(byrs_tot, byrs_tot);

  // Year total span from the first birth year to the last sample year e.g., samples years 2015-2018, birth years 2010,2017
  // The cumulative survival is computed for years 2010-2018. Remember 0 index
  // The index for the sampling years is yrs_tot - syrs_tot  to syrs_tot. Check yrs_tot = (2018 - 2010) = 8
  // The index for the birth years is yrs_tot - syrs_tot  to syrs_tot. Check yrs_tot = (2018 - 2010) = 8

  // ===============
  // Female dynamics
  // ===============
  N_f[0] = exp(N0_f);
  for (year = 1; year < yrs_tot; year++)
  {
    N_f[year] = N_f[0] * exp(R1 * year);
  }
    
  // =============
  // Male dynamics 
  // =============
  // N_m[0] = N0_m;
  // for (year = 1; year < yrs_tot; year++)
  // {
  //   N_m[year] = N0_m * exp(R2 * year);
  // }

  // ======================
  // CKMR POP probabilities
  // ======================
   int Yad, Bju, rel_Yad;
   for (Yad = 0; Yad < syrs_tot; Yad++)
   {
     rel_Yad = Yad + start_yrs_diff;
     for (Bju = 0; Bju < byrs_tot; Bju++)
     {
       // Female survival
       Type cumul_surv_f1 = exp(-Z1 * (rel_Yad - Bju));
       Type cumul_surv_f2 = exp(-Z1 * (Bju - rel_Yad));
       // Male survival
//       Type cumul_surv_m1 = exp(-Z2 * (rel_Yad - Bju));
//       Type cumul_surv_m2 = exp(-Z2 * (Bju - rel_Yad));
      
       std::cout << "Proposed index " << rel_Yad << std::endl;
       std::cout << "Indices " << Yad  << " " << Bju << std::endl;
       if (Bju < rel_Yad)
       {
         std::cout << "POP prob 1" << std::endl;
         mop_prob(Bju, Yad) = cumul_surv_f1 / N_f[rel_Yad];
//         fop_prob(Bju, Yad) = cumul_surv_m1 / N_m[rel_Yad];
       } else if (!lethal && (Bju > rel_Yad))
       {
         std::cout << "POP prob 2" << std::endl;
         mop_prob(Bju, Yad) = cumul_surv_f2 / N_f[Bju];
//         fop_prob(Bju, Yad) = cumul_surv_m2 / N_m[Bju];
       }
     }
  }

   std::cout << "MOP probabilities" << std::endl;
   for (year = 0; year < byrs_tot; year++)
   {
     std::cout << mop_prob.row(year) << std::endl;
   }

  // ======================
  // CKMR HSP probabilities
  // ======================
  for (birth_year_1 = 0; birth_year_1 < byrs_tot; birth_year_1++)
  {
    for (birth_year_2 = 0; birth_year_2 < byrs_tot; birth_year_2++)
    {
      if ((birth_year_1 <= birth_year_2))
      {
		// Female survival
        Type cumul_surv_f = exp(-Z1 * (birth_year_2 - birth_year_1));
        // Male survival
        // Type cumul_surv_m = exp(-Z2 * (birth_year_2 - birth_year_1));
		
        mhsp_prob(birth_year_1, birth_year_2) = cumul_surv_f / N_f[birth_year_2];
        // fhsp_prob(birth_year_1, birth_year_2) = cumul_surv_m / N_m[birth_year_2];
      }
    }
  }

  mhsp_prob = mhsp_prob * fudge_factor;
  //std::cout << "MOP probabilities" << std::endl;
  //std::cout << mop_prob.row(0) << std::endl;
  //std::cout << "FOP probabilities" << std::endl;
  //std::cout << fop_prob.row(0) << std::endl;
  std::cout << "MHSP probabilities" << std::endl;
  std::cout << mhsp_prob.row(0) << std::endl;
  //std::cout << "FHSP probabilities" << std::endl;
  //std::cout << fhsp_prob.row(0) << std::endl;
	  
  
  // ==================
  // Log-likelihood POP
  // ==================
//  std::cout << "Starting POP llk" << std::endl;
//  std::cout << "Sampling years total"  << syrs_tot << std::endl;
//  std::cout << "Birth years total" << byrs_tot<< std::endl;
  ciff_llk = 0.0;
   for (Yad = 0; Yad < syrs_tot; Yad++)
   {
	 rel_Yad = Yad + start_yrs_diff;
     for (Bju = 0; Bju < byrs_tot; Bju++)
     {
	   if (Bju < rel_Yad || (!lethal && Bju > rel_Yad))
	   {
	 	  ciff_llk -= dbinom(kin_mop(Bju, Yad),
	 						 comb_mop(Bju, Yad),
	 						 mop_prob(Bju, Yad), true);
          exp_mops(Bju, Yad) = comb_mop(Bju, Yad) * mop_prob(Bju, Yad);
	 	  std::cout << "MOP llk " << ciff_llk << std::endl;
	 	  std::cout << "year birth year " << Yad << " " << Bju <<
	 				 " kin " << kin_mop(Bju, Yad) <<
	 				 " combs " << comb_mop(Bju, Yad) <<
	 				 " prob " << mop_prob(Bju, Yad) << std::endl;
//	 	  ciff_llk -= dbinom(kin_fop(Bju, Yad),
//	 						 comb_fop(Bju, Yad),
//	 						 fop_prob(Bju, Yad), true);
	 //      std::cout << "FOP llk " << ciff_llk << std::endl;
	 //      std::cout << "year birth year " << Yad << " " << Bju << std::endl;
	   }
     }
   }


  // ==================
  // Log-likelihood HSP
  // ==================
  std::cout << "Starting HSP llk" << std::endl;
  for (birth_year_1 = 0; birth_year_1 < byrs_tot; birth_year_1++)
  {
    for (birth_year_2 = 0; birth_year_2 < byrs_tot; birth_year_2++)
    {
      if (birth_year_1 <= birth_year_2)
      {
		std::cout << "HSP llk " << ciff_llk << std::endl;
		std::cout << "by1 by2 " << birth_year_1 << " " << birth_year_2 << std::endl;
        ciff_llk -= dbinom(kin_mhsp(birth_year_1, birth_year_2),
                           comb_mhsp(birth_year_1, birth_year_2),
                           mhsp_prob(birth_year_1, birth_year_2), true);
        exp_mhsps(birth_year_1, birth_year_2) = comb_mhsp(birth_year_1, birth_year_2) * mhsp_prob(birth_year_1, birth_year_2);
        // ciff_llk -= dbinom(kin_fhsp(birth_year_1, birth_year_2),
        //                    comb_fhsp(birth_year_1, birth_year_2),
        //                    fhsp_prob(birth_year_1, birth_year_2), true);
      }
    }
  }

  REPORT(mop_prob);
  REPORT(mhsp_prob);
  REPORT(exp_mops);
  REPORT(exp_mhsps);
  ADREPORT(N_f);
  ADREPORT(log(N_f));
  return ciff_llk;
}
