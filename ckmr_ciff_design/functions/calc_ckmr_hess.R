# Assumes the function environment has been set to com.kin

calc_ckmr_hess <- function(pars)
{
  # pars is vector of starting parameters (Nf_0, Nm_0, Z1, Z2, R1, R2)
  N0_f <- pars[1]
  N0_m <- pars[2]
  Z1   <- pars[3]
  Z2   <- pars[4]
  R1   <- pars[5]
  R2   <- pars[6]
  
  # Set up comparison arrays
  dim_bju <- dim(n_comp_MOP_noage)[1]
  dim_yad <- dim(n_comp_MOP_noage)[2]
  pop_comp_array <- array(0, dim= c(dim_bju, dim_yad, 2))
  
  pop_comp_array[,,1] <- as.array(n_comp_MOP_noage)       # MOP comparisons
  pop_comp_array[,,2] <- as.array(n_comp_FOP_noage)       # FOP comparisons
  hsp_comp_array      <- as.array(n_comp_HSP)             # HSP comparisons (Both MHSP and FHSP)
  
  hess_temp <- matrix(0, nrow = length(pars), ncol = length(pars))
  # For all birth years
  for (i in seq(1, dim_bju))
  {
    # Consider all adult sampling years
    for (j in seq(1, dim_yad))
    {
      # Look at both MOPs (k=1) and FOPs (k=2)
      for (k in seq(1, 2))
      {
        data <- list(symin   = min(Yad_range), symax = max(Yad_range), 
                     bymin   = min(Bju_range), bymax = max(Bju_range), 
                     lethal  = lethal, pop_bju = i,   pop_yad = j,
                     mode = "PROB",    setting = k)
        ciff_prob <- MakeADFun(data, 
                               parameters = list(N0_f=N0_f, N0_m=N0_m, Z1=Z1, Z2=Z2, R1=R1, R2=R2), silent=TRUE)
        prob_grad_kin_zz <- ciff_prob$gr(pars)
        h_kin_zz         <- 4 * t(prob_grad_kin_zz) %*% prob_grad_kin_zz
        hess_temp        <- hess_temp +  pop_comp_array[i, j, k] * h_kin_zz
      }
    }
    
    # Consider all birth years (for HS)
    for (j in seq(1, dim_bju)) 
    {
      # Look at both MHSPs (k=3) and FHSPs (k=4)
      for (k in seq(3, 4))
      {
        data <- list(symin   = min(Yad_range), symax = max(Yad_range),
                     bymin   = min(Bju_range), bymax = max(Bju_range),
                     lethal  = lethal, pop_bju = i, pop_yad = j, # Here Yad is actually Bju of older HS
                     mode = "PROB",    setting = k)
        ciff_prob <- MakeADFun(data,
                               parameters = list(N0_f=N0_f, N0_m=N0_m, Z1=Z1, Z2=Z2, R1=R1, R2=R2), silent=TRUE)
        prob_grad_kin_zz <- ciff_prob$gr(pars)
        h_kin_zz         <- 4 * t(prob_grad_kin_zz) %*% prob_grad_kin_zz
        hess_temp        <- hess_temp +  hsp_comp_array[i, j] * h_kin_zz
      }
    }
  }
  
  # Return hessian / sum of all combinations
  return(hess_temp)
}