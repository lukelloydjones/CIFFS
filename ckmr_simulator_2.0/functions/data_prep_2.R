data_prep_2 <- function(sam.no.py, years.samp, scen, repi, lethal, two_or_none)
{
  # PARS
  #  sam.no.py  - number of individuals to sample from each year
  #  years.samp - number of years to sample
  #  scen       - what simulation scenario to sample. Defined by mortality parameter
  #  repi - what simulation rep are we doing
  #  lethal - 1 => Lethal sampling, 0 => mark-recapture
  # RETURNS
  #  saves the environment to be passed to log likelihood for fitting
  #      com.kin contains
  #       -
  # ---------------------------------------
  # Make a matrix of size n.x.n and compare
  # each individual
  # For each comparison record the relationship if it exists
  # ---------------------------------------
  source("functions/prep_comb_kin_2.R")
  print(paste0("Preping kin and combinations for population sample ", repi))
  # ----------------------
  # Making of combinations
  # ----------------------
  load(paste0("found_kin", two_or_none ,"/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
              years.samp, "_samp_no_py_", sam.no.py, ".Rdata"))
  pop.samp.i  <- read.table(paste0("sampled_populations", two_or_none, "/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
                            years.samp, "_samp_no_py_", sam.no.py, ".txt"), 
                            header = T)
  com.kin <- prep_comb_kin_2(kin, pop.samp.i, lethal)
  save(com.kin, file = paste0("preped_data", two_or_none, "/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
                              years.samp, "_samp_no_py_", sam.no.py, ".Rdata"))
}