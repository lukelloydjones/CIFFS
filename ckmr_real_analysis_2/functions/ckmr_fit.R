ckmr_fit <- function(sam.no.py, years.samp, scen, repi)
{
  load(paste0("preped_data/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
         years.samp, "_samp_no_py_", sam.no.py, ".Rdata"))
  print(class(com.kin))
  print(ls(com.kin))
  environment(ciff_ckmr_loglik) <- com.kin
  # -----------
  # Numeric fit
  # -----------
  fit1 = nlminb(c(1000, 1000, 0.4, 0.4, 1.1, 1.1), ciff_ckmr_loglik,
                lower = c(1, 1, -10, -10, -10 , -11),
                upper = c(10000, 10000, 10, 10, 200, 200))
  print(ls(com.kin))
  # print(fit1$par)
  # N <- com.kin$N
  # pars <- fit1$par
  # print("Abundance through time")
  # print(com.kin$N)
  # print("No. POPs")
  # print(sum(com.kin$n_MOP_noage) + sum(com.kin$n_FOP_noag))
  # print("No. HSPs")
  # print(sum(com.kin$n_mhsps) + sum(com.kin$n_phsps))
  # return(list(N, pars))
}