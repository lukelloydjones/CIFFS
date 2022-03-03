# =============================================================================
# Outline of requirements for CKMR design basics
# =============================================================================
args = commandArgs(trailingOnly=TRUE)
max_samp_size <- as.numeric(args[1])
juv_max_samp  <- as.numeric(args[2])
sad_max_samp  <- as.numeric(args[3])
adu_max_samp  <- as.numeric(args[4])
res_type      <- as.numeric(args[5])


# Then we need a TMB function to compute the first order partial derivatives 
# with repsect to all paramters in the model for all covariate types ewwwww
suppressPackageStartupMessages({
  library(atease)   
  library(mvbutils)
  library(dplyr)
  library(offarray)
  #library(ggplot2)
  #library(reshape)
  #library(RColorBrewer)
  library(TMB)
  library(rgenoud)
})
make_directories <- function(...) {
  dirs <- unlist(list(...))
  stopifnot(length(dirs) > 0, is.character(dirs))
  for(d in dirs) {
    if(!dir.exists(d)) dir.create(d)
  }
  invisible(d)
}

# Compile the AD components
compile("ciff_cpp/ciff_prob_deriv.cpp")
dyn.load(dynlib("ciff_cpp/ciff_prob_deriv")) 

source("functions/function_list.R")

# ---------------------------------
# Partition the data in myriad ways
# ---------------------------------
make_directories(file.path("preped_data"))
#ex.dir <- "/Users/llo080/Dropbox/CSIRO/summer_scholarship/smp_summer_isabel_peters/ckmr_simulator_2.0"
#ex.dir <- "/Users/pet247/Dropbox/smp_summer_isabel_peters/ckmr_simulator_2.0"
years.samp <- 5    # How many years do you want to sample

# ---------------------------------
# Now make a function
# ---------------------------------

ckmr_design_func <- function(theta, type, max_total_samps = 0)
{
  # print(prob_vec)
  # print(type)
  # PARS
  #  theta[1:years.samp] - vector of proportion of individuals to sample in each of 1:years.samp years
  #  theta[years.samp+1] - prop_ju_py - proportion of juveniles in each year
  #  theta[years.samp+2] - prop_sa_py - proportion of aub-adults in each year
  #  theta[years.samp+3] - prop_ad_py - proportion of adulds in each year
  #  type                - what type of objective function to return 1 -min CV, 2 - sum CV over years, 3 the CV vector for plotting
  # RETURNS
  #  A variety of CV related objective function values depends on type 
  #  abund_cv            - cv of the abundance at each year
  
  tryCatch({
  
  samp.no.py <- matrix(floor(theta), nrow = years.samp, ncol = 3, byrow = T)
  if (max_total_samps != 0 && sum(samp.no.py) > max_total_samps)
  {
    stop("Max total samps exceeded")
  }
  dimnames(samp.no.py) <- list(c(2014:2018), c('ju', 'subad', 'ad'))
  
  #rep(c(90, 50), 5)  # How many individuals do you want to sample per year
  lethal        <- 0    # 1 => Lethal sampling, 0 => mark-recapture
  scen          <- 0.3  # Of the scenarios simulated, which one do you want to analayse?
  
  f_prop = 0.45
  com.kin <- prep_comb_kin_2_var_props(samp.no.py, years.samp, f_prop, lethal)
  
  # ----------------------------------------------------------
  # Calculate the CKMR approximate to the LLGK hessian
  # ----------------------------------------------------------
  # Need a set of starting parameters Nf_0, Nm_0, Z1, Z2, R1, R2
  pars <- c(800,  874,  0.32,  0.29,  -0.045,  -0.050) # Hypothesised truth
  N0_f <- pars[1]
  N0_m <- pars[2]
  Z1   <- pars[3]
  Z2   <- pars[4]
  R1   <- pars[5]
  R2   <- pars[6]
  
  environment(calc_ckmr_hess) <- com.kin
  ckmr_hess <- calc_ckmr_hess(pars)
  
  # ----------------------------------------------------------
  # Calculate the CVs for the abundance estimates through time
  # ----------------------------------------------------------
  
  abund_cv <- array(0, length(com.kin$years))
  names(abund_cv) <- as.character(com.kin$years)
  for (year in com.kin$years) {
    data  <- list(symin    = min(com.kin$Yad_range), symax = max(com.kin$Yad_range),
                  bymin    = min(com.kin$Bju_range), bymax = max(com.kin$Bju_range),
                  lethal   = com.kin$lethal,
                  pop_bju = 1, pop_yad = 1,
                  mode = "ABUND", setting = year)
    ciff_abund_grad <- MakeADFun(data,
                                 parameters = list(N0_f=N0_f, N0_m=N0_m, Z1=Z1, Z2=Z2, R1=R1, R2=R2), silent=TRUE)
    abund_grad_year <- ciff_abund_grad$gr(pars)
    N <- ciff_abund_grad$report()$N
    names(N) <- as.character(com.kin$years)
    var             <- abund_grad_year %*% solve(ckmr_hess) %*% t(abund_grad_year)
    cv              <- sqrt(var) / N[as.character(year)]
    abund_cv[as.character(year)] = cv
  }
  # Compute expected kin
  hsps <- ciff_abund_grad$report()$mhsp_prob * com.kin$n_comp_HSP + ciff_abund_grad$report()$fhsp_prob * com.kin$n_comp_HSP
  pops <- ciff_abund_grad$report()$mop_prob * com.kin$n_comp_MOP_noage + ciff_abund_grad$report()$fop_prob * com.kin$n_comp_FOP_noage
  
  # Return different objective functions
  print("****************************************")
  print("Samps py")
  print(samp.no.py)
  print("MIN CV, SUM OF CVs, CV vector")
  print(c(min(abund_cv), sum(abund_cv)))
  print(abund_cv)
  print("Expected number of POPs and HSPs")
  print(c(sum(pops), sum(hsps)))
  print("****************************************")
  
  # ---------------- 
  # Draw a nice plot
  # ----------------
  # plot_abund_cv(com.kin$years, N, abund_cv, lethal=0, 100, 3000)
  
  # What type of objective function
  if (type == 1)
  {
    return(min(abund_cv))
  } else if (type == 2)
  {
    return(sum(abund_cv))
  } else {
    return(abund_cv)
  }
  
  
  #return(sum(abund_cv) + diff)
  #return(list(min(abund_cv), sum(abund_cv), abund_cv)) # Minimum over CV vector, Sum of CV vector
  
  }, error=function(e) {
    print(paste0("ERROR - ", geterrmessage()))
    return(100)
  })
}

# ---------------------
# Now optimise
# ---------------------
d <- matrix(rep(c(0, juv_max_samp, 0, sad_max_samp, 0, adu_max_samp), 5), nrow = 15, ncol = 2, byrow = T)
out_genoud <- genoud(ckmr_design_func, nvars=15, pop.size=80, max.generations=80, wait.generations=10,
                     hard.generation.limit=F, starting.values=rep(10, 15),
                     Domains=d,solution.tolerance=1e-4,
                     boundary.enforcement=1, 
                     data.type.int=T, 
                     type = res_type, max_total_samps = max_samp_size)

out_pars <- matrix(out_genoud$par, nrow = years.samp, ncol = 3, byrow = T)
dimnames(out_pars) <- list(c(2014:2018), c('ju', 'subad', 'ad'))


out_file_no = args[6]
write.table(out_pars, paste("out/design", max_samp_size, res_type, out_file_no, "pp
ar.txt", sep = "_"), col.names = F,
            row.names = F, quote = F, sep = "\t")
write.table(out_genoud$value, paste("out/design", max_samp_size, res_type, out_file_no,  
"val.txt", sep = "_"), col.names = F,
            row.names = F, quote = F, sep = "\t")




