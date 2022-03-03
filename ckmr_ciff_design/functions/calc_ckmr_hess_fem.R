# Assumes the function environment has been set to com.kin
# Testing 
#n_comp_MOP_noage <- com.kin$n_comp_MOP_noage
#n_comp_HSP <- com.kin$n_comp_HSP
#Yad_range  <- com.kin$Yad_range
#Bju_range  <- com.kin$Bju_range
calc_ckmr_hess_fem <- function(pars)
{
  # pars is vector of starting parameters (Nf_0, Nm_0, Z1, Z2, R1, R2)
  N0_f <- pars[1]
  Z1   <- pars[3]
  R1   <- pars[5]
  
  # Set up comparison arrays
  dim_bju <- dim(n_comp_MOP_noage)[1]
  dim_yad <- dim(n_comp_MOP_noage)[2]
  pop_comp_array <- array(0, dim= c(dim_bju, dim_yad, 2))
  
  pop_comp_array[,,1] <- as.array(n_comp_MOP_noage)       # MOP comparison
  hsp_comp_array      <- as.array(n_comp_HSP)             # HSP comparisons (Both MHSP and FHSP)
  
  hess_temp <- matrix(0, nrow = length(pars), ncol = length(pars))
  # For all birth years
  for (i in seq(1, dim_bju))
  {
    # Consider all adult sampling years
    for (j in seq(1, dim_yad))
    {
      # Look at both MOPs (k=1) and FOPs (k=2)
        data <- list(symin   = min(Yad_range), symax = max(Yad_range), 
                     bymin   = min(Bju_range), bymax = max(Bju_range), 
                     lethal  = lethal, pop_bju = i,   pop_yad = j,
                     mode = "PROB",    setting = 1)
        ciff_prob <- MakeADFun(data, 
                               parameters = list(N0_f=N0_f, Z1=Z1, R1=R1), silent=TRUE, DLL = "ciff_prob_deriv_fem")
        prob_grad_kin_zz <- ciff_prob$gr(pars)
        h_kin_zz         <- 4 * t(prob_grad_kin_zz) %*% prob_grad_kin_zz
        hess_temp        <- hess_temp +  pop_comp_array[i, j, 1] * h_kin_zz
    }
    
    # Consider all birth years (for HS)
    for (j in seq(1, dim_bju)) 
    {
      # Look at both MHSPs (k=3) and FHSPs (k=4)
        data <- list(symin   = min(Yad_range), symax = max(Yad_range),
                     bymin   = min(Bju_range), bymax = max(Bju_range),
                     lethal  = lethal, pop_bju = i, pop_yad = j, # Here Yad is actually Bju of older HS
                     mode = "PROB",    setting = 2)
        ciff_prob <- MakeADFun(data,
                               parameters = list(N0_f=N0_f, Z1=Z1, R1=R1), silent=TRUE, DLL = "ciff_prob_deriv_fem")
        prob_grad_kin_zz <- ciff_prob$gr(pars)
        h_kin_zz         <- 4 * t(prob_grad_kin_zz) %*% prob_grad_kin_zz
        hess_temp        <- hess_temp +  hsp_comp_array[i, j] * h_kin_zz
    }
  }
  
  # Return hessian / sum of all combinations
  return(hess_temp)
}