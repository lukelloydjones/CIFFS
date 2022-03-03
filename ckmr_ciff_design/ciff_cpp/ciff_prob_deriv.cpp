#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  //DATA_MATRIX(comb_mop);                // Data vector transmitted from R
  //DATA_MATRIX(comb_fop);                // Data vector transmitted from R
  //DATA_MATRIX(comb_mhsp);               // Data vector transmitted from R
  //DATA_MATRIX(comb_fhsp);               // Data vector transmitted from R
    
  DATA_INTEGER(symin);                  // Data vector transmitted from R
  DATA_INTEGER(symax);                  // Data vector transmitted from R
  DATA_INTEGER(bymin);                  // Data vector transmitted from R
  DATA_INTEGER(bymax);                  // Data vector transmitted from R
  DATA_INTEGER(lethal);                 // Data vector transmitted from R
  
  DATA_INTEGER(pop_bju);                // Data vector transmitted from R
  DATA_INTEGER(pop_yad);                // Data vector transmitted from R
  
  DATA_STRING(mode);                    // Mode transmitted from R. Either PROB or ABUND
  DATA_INTEGER(setting);                // Setting transmitted from R. 
  // In PROB mode this is kin_type, in ABUND mode this is abund_year
  
  // Parameter values transmitted from R
  PARAMETER(N0_f);                      // - Initial number of females
  PARAMETER(N0_m);                      // - Initial number of males
  PARAMETER(Z1);                        // - Female mortality
  PARAMETER(Z2);                        // - Male mortality
  PARAMETER(R1);                        // - Female recruitment
  PARAMETER(R2);                        // - Male recruitment


  // Initial declarations
  int yrs_tot, byrs_tot, syrs_tot, year, birth_year_1, birth_year_2, start_yrs_diff;
  yrs_tot  = symax - bymin + 1;
  syrs_tot = symax -  symin + 1;
  byrs_tot = bymax - bymin + 1;
  start_yrs_diff = symin - bymin;

  // std::cout << "Years total " << yrs_tot << std::endl;
  // std::cout << "syrs_tot " << syrs_tot << std::endl;
  // std::cout << "byrs_tot" << byrs_tot << std::endl;
    
  vector<Type> N_f(yrs_tot);
  vector<Type> N_m(yrs_tot);
  vector<Type> N(yrs_tot);
  
  matrix<Type> mop_prob(byrs_tot, syrs_tot);
  matrix<Type> fop_prob(byrs_tot, syrs_tot);
  matrix<Type> mhsp_prob(byrs_tot, byrs_tot);
  matrix<Type> fhsp_prob(byrs_tot, byrs_tot);

  // Year total span from the first birth year to the last sample year 
  // e.g., samples years 2015-2018, birth years 2010-2017
  // The cumulative survival is computed for years 2010-2018. Remember 0 index
  // The index for the sampling years is yrs_tot - syrs_tot  to syrs_tot.
  // The index for the birth years is yrs_tot - byrs_tot  to byrs_tot. 
  // Check yrs_tot = (2018 - 2010) = 8

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
  N_m[0] = exp(N0_m);
  for (year = 1; year < yrs_tot; year++)
  {
    N_m[year] = N_m[0] * exp(R2 * year);
  }
  
  N = N_f + N_m;
  REPORT(N);
  ADREPORT(log(N));
  
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
      Type cumul_surv_m1 = exp(-Z2 * (rel_Yad - Bju));
      Type cumul_surv_m2 = exp(-Z2 * (Bju - rel_Yad));
      
      //std::cout << "Proposed index " << rel_Yad << std::endl;
      //std::cout << "Indices " << Yad  << " " << Bju << std::endl;
      if (Bju < rel_Yad)
      {
        //std::cout << "POP prob 1" << std::endl;
        mop_prob(Bju, Yad) = cumul_surv_f1 / N_f[rel_Yad];
        fop_prob(Bju, Yad) = cumul_surv_m1 / N_m[rel_Yad];
      } else if (!lethal && (Bju > rel_Yad))
      {
        //std::cout << "POP prob 2" << std::endl;
        mop_prob(Bju, Yad) = cumul_surv_f2 / N_f[Bju];
        fop_prob(Bju, Yad) = cumul_surv_m2 / N_m[Bju];
      }
    }
  }

  // std::cout << "MOP probabilities" << std::endl;
  // for (year = 0; year < byrs_tot; year++)
  // {
    // std::cout << mop_prob.row(year) << std::endl;
  // }

  // ======================
  // CKMR HSP probabilities
  // ======================
  for (birth_year_1 = 0; birth_year_1 < byrs_tot; birth_year_1++)
  {
    for (birth_year_2 = 0; birth_year_2 < byrs_tot; birth_year_2++)
    {
      if ((birth_year_1 < birth_year_2))
      {
	 	// Female survival
        Type cumul_surv_f = exp(-Z1 * (birth_year_2 - birth_year_1));
        // Male survival
        Type cumul_surv_m = exp(-Z2 * (birth_year_2 - birth_year_1));

        mhsp_prob(birth_year_1, birth_year_2) = cumul_surv_f / N_f[birth_year_2];
        fhsp_prob(birth_year_1, birth_year_2) = cumul_surv_m / N_m[birth_year_2];
      }
    }
  }

  REPORT(mhsp_prob);
  REPORT(fhsp_prob);
  REPORT(mop_prob);
  REPORT(fop_prob);
   
  if (mode == "PROB") 
  {
	int kin_type = setting;
	if (kin_type == 1)
	{
	  return sqrt(mop_prob(pop_bju-1, pop_yad-1));
	} else if  (kin_type == 2)
	{
	  return sqrt(fop_prob(pop_bju-1, pop_yad-1));
	} else if  (kin_type == 3)
	{
	  return sqrt(mhsp_prob(pop_bju-1, pop_yad-1));
	} else if  (kin_type == 4)
	{
	  return sqrt(fhsp_prob(pop_bju-1, pop_yad-1));
	} else {
		std::cout << "PROB mode: Invalid kin type provided" << std::endl;
	}
	
  } else if (mode == "ABUND") 
  {
	int abund_year = setting - bymin;
    if ((0 <= abund_year) && (abund_year < yrs_tot))
	{
	  return (N_f[abund_year] + N_m[abund_year]);
	} else {
		std::cout << "ABUND mode: Invalid abund year provided" << std::endl;
	}
	
  } else if (mode == "LOGABUND")
  {
    int abund_year = setting - bymin;
    if ((0 <= abund_year) && (abund_year < yrs_tot))
    {
      return (log(N_f[abund_year] + N_m[abund_year]));
    } else {
        std::cout << "LOGABUND mode: Invalid abund year provided" << std::endl;
    }
    
  } else {
	std::cout << "Invalid mode provided" << std::endl;
  }
}
