# =============================================================================
# Outline of requirements for CKMR design basics
# =============================================================================
# Then we need a TMB function to compute the first order partial derivatives 
# with repsect to all paramters in the model for all covariate types ewwwww
suppressPackageStartupMessages({
  library(atease)   
  library(mvbutils)
  library(dplyr)
  library(offarray)
  library(ggplot2)
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

ckmr_design_func <- function(theta, type, max_total_samps = 0, draw_plot = 0, min = 0, max = 5000)
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
  dimnames(samp.no.py) <- list(c(2015:(2014+years.samp)), c('ju', 'subad', 'ad'))
  
  #rep(c(90, 50), 5)  # How many individuals do you want to sample per year
  lethal        <- 0    # 1 => Lethal sampling, 0 => mark-recapture
  scen          <- 0.3  # Of the scenarios simulated, which one do you want to analayse?
  
  f_prop = 0.35
  com.kin <- prep_comb_kin_2_var_props(samp.no.py, years.samp, f_prop, lethal)
  
  # ----------------------------------------------------------
  # Calculate the CKMR approximate to the LLGK hessian
  # ----------------------------------------------------------
  # Need a set of starting parameters Nf_0, Nm_0, Z1, Z2, R1, R2
  pars <- c(log(800), log(874), 0.32, 0.29, -0.045, -0.050) # Hypothesised truth
  pars <- c(log(800), 0, 0.32, 0, -0.045, 0) # Hypothesised truth
  #800, 874, 0.32, 0.29, -0.045, -0.050
  #600, 600, 1.79, 1.79, -1,-1
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
  std.err  <- array(0, length(com.kin$years))
  names(abund_cv) <- as.character(com.kin$years)
  names(std.err)  <- as.character(com.kin$years)
  for (year in com.kin$years) {
    data  <- list(symin    = min(com.kin$Yad_range), symax = max(com.kin$Yad_range),
                  bymin    = min(com.kin$Bju_range), bymax = max(com.kin$Bju_range),
                  lethal   = com.kin$lethal,
                  pop_bju = 1, pop_yad = 1,
                  mode = "LOGABUND", setting = year)
    ciff_abund_grad <- MakeADFun(data,
                                 parameters = list(N0_f=N0_f, N0_m=N0_m, Z1=Z1, Z2=Z2, R1=R1, R2=R2), silent=TRUE)
    abund_grad_year <- ciff_abund_grad$gr(pars)
    rep  <- sdreport(ciff_abund_grad)
    logN <- summary(rep)[grep("log.", rownames(summary(rep))), 1]
    names(logN) <- as.character(com.kin$years)
    var      <- abund_grad_year %*% solve(ckmr_hess) %*% t(abund_grad_year)
    std.err[as.character(year)] <- sqrt(var)
    cv       <- sqrt(exp(var) - 1)
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
  if (draw_plot == 1) {
    #plot_abund_cv(com.kin$years, N, abund_cv, lethal=0)
    #plot_abund_cv(com.kin$years, N, abund_cv, lethal=0, 100, 3000)
    #plot_abund_cv(com.kin$years, N, abund_cv, lethal=0, 0, 3100)
    rep <- sdreport(ciff_abund_grad)
    N.res <- as.data.frame(summary(rep)[grep("log.", rownames(summary(rep))), 1])
    # Build the results into a item to pass to the plot function
    plot_abund_cv(com.kin$years, N.res, abund_cv, std.err, lethal = 0, min, max)
  }
  
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

# ckmr_design_func(c(rep(100,5), 300, 300), 2)
# ckmr_design_func(c(rep(50,5), 532, 30), 1)
# 
# ckmr_design_func(c(9, 76, 0, 19, 61, 0, 31, 72, 0, 42, 62, 4, 76, 39, 9), 1)
# ckmr_design_func(c(4, 89, 1, 30, 65, 0, 58, 25, 2, 31, 66, 5, 37, 87, 0), 1)
# ckmr_design_func(c(9, 52, 6, 20, 41, 0, 42, 4,  0, 69, 19, 0, 91, 90, 57), 1)
# 
# rep(c(50, 50, 5),5)
# ckmr_design_func(rep(c(50, 50, 5),5), 1)
# 
# 
# nlminb(theta, ckmr_design_func, gradient = NULL, hessian = NULL, age.s1, 2,
#        control = list(step.min = 100),
#        lower = c(rep(0.01, 10), c(0.01, 0.01, 0.01)),
#        upper = c(rep(1, 10),    c(1, 1, 1)))



# -------------------------------------------------------
# Some experimentation
# -------------------------------------------------------

ckmr_design_func(rep(c(20, 20, 60), 5), 1, draw_plot = 1)


ckmr_design_func(rep(c(30, 30, 15), 5), 1, draw_plot = 1)
ckmr_design_func(rep(c(40, 40, 20), 5), 1, draw_plot = 1)
ckmr_design_func(rep(c(50, 50, 25), 5), 1, draw_plot = 1)

ckmr_design_func(rep(c(60, 60, 30), 5), 1, draw_plot = 1)

ckmr_design_func(rep(c(25, 25, 10), 5), 1, draw_plot = 1)
ckmr_design_func(rep(c(50, 50, 20), 5), 1, draw_plot = 1)



###########################################################
# All 728 inds
ckmr_design_func(c(68, 10, 27, 
                   115, 32, 56, 
                   77, 26, 55, 
                   79, 2, 31, 
                   77, 9, 64), 1, draw_plot = 1, min = 500, max = 4000)
# 7 plates # Remove 56 adults. 672 individuals
ckmr_design_func(c(68,  10, 0, 
                   115, 32, 27, 
                   77,  26, 55, 
                   79,  2, 31, 
                   77,  9, 64), 1, draw_plot = 1, min = 500, max = 4000)
# 6 plates # Remove 96 adults. 576 individuals
ckmr_design_func(c(68,  10, 0, 
                   115, 32, 0, 
                   77,  26, 20, 
                   79,  2, 20, 
                   77,  9, 41), 1, draw_plot = 1, min = 500, max = 4000)
# 5 plates, 480 samples
ckmr_design_func(c(63, 10,  0, 
                   105, 32, 0, 
                   77, 26, 0, 
                   79, 2, 0, 
                   77, 9, 0), 1, draw_plot = 1, min = 500, max = 4000)
# 4 plates, 384 samples
ckmr_design_func(c(43,  10, 0, 
                   69,  32, 0, 
                   57,  26, 0, 
                   59,  2,  0, 
                   77,  9,  0), 1, draw_plot = 1, min = 500, max = 4000)
# 3 plates, 288 samples
ckmr_design_func(c(20,  10, 0, 
                   30,  32, 0, 
                   37,  26, 0, 
                   49,  2,  0, 
                   73,  9,  0), 1, draw_plot = 1, min = 500, max = 5000)

ckmr_design_func(c(32,  0,  6, 55,  4,  9,  4, 15,  7, 33,  6,  3, 52, 28, 34), 1, draw_plot = 1)


d <- matrix(c(rep(0, 15), c(68, 10, 27, 
                            115, 32, 56, 
                            77, 26, 55, 
                            79, 2, 31, 
                            77, 9, 64)), nrow = 15, ncol = 2)
out_genoud <- genoud(ckmr_design_func, nvars=15, pop.size=80, max.generations=80, wait.generations=10,
                     hard.generation.limit=F,
                     Domains=d,solution.tolerance=1e-4,
                     boundary.enforcement=1, 
                     data.type.int=T, 
                     starting.values = c(43,  10, 0, 
                                         69,  32, 0, 
                                         57,  26, 0, 
                                         59,  2,  0, 
                                         77,  9,  0),
                     type = 1, max_total_samps = 384)

out_pars <- matrix(out_genoud$par, nrow = years.samp, ncol = 3, byrow = T)
dimnames(out_pars) <- list(c(2015:(2014+years.samp)), c('ju', 'subad', 'ad'))
out_pars


###########################################################

# 45, 0, 18, 65, 8, 49, 8, 20, 15, 79, 17, 17, 78, 34, 46

ckmr_design_func(c(45, 0, 18, 65, 8, 49, 8, 20, 15), 1, draw_plot = 1)

ckmr_design_func(c(45, 0, 18, 65, 8, 49, 8, 20, 15, 79, 17, 17, 78, 34, 27), 1, draw_plot = 1)
ckmr_design_func(c(45, 0,  2, 65, 8,  5, 8, 20,  2, 79, 17,  5, 78, 34, 16), 1, draw_plot = 1)
ckmr_design_func(c(19, 0,  0, 25, 8,  0, 8, 20,  0, 79, 17,  0, 78, 34,  0), 1, draw_plot = 1)

ckmr_design_func(c(45, 0, 12, 61, 8, 48, 8, 20, 15, 79, 16, 17, 76, 34, 41), 1, draw_plot = 1)
ckmr_design_func(c(45, 0,  0, 50, 6,  5, 8, 20,  2, 79, 17,  5, 78, 34, 35), 1, draw_plot = 1)
ckmr_design_func(c(36, 0,  0, 20, 8,  0, 8, 19,  0, 79, 14,  3, 74, 24,  3), 1, draw_plot = 1)
# Much improved

# d <- matrix(rep(c(0, 100, 0, 100, 0, 50), 5), nrow = 15, ncol = 2, byrow = T)
d <- matrix(c(rep(0, 15), 45, 0, 18, 65, 8, 49, 8, 20, 15, 79, 17, 17, 78, 34, 46), nrow = 15, ncol = 2)
out_genoud <- genoud(ckmr_design_func, nvars=15, pop.size=80, max.generations=80, wait.generations=10,
                     hard.generation.limit=F,
                     Domains=d,solution.tolerance=1e-4,
                     boundary.enforcement=1, 
                     data.type.int=T, 
                     starting.values = c(45, 0, 12, 61, 8, 48, 8, 20, 15, 79, 16, 17, 76, 34, 41),
                     type = 1, max_total_samps = 480)

out_pars <- matrix(out_genoud$par, nrow = years.samp, ncol = 3, byrow = T)
dimnames(out_pars) <- list(c(2015:(2014+years.samp)), c('ju', 'subad', 'ad'))

diff_pars <- matrix(c(45, 0, 18, 65, 8, 49, 8, 20, 15, 79, 17, 17, 78, 34, 46) - out_genoud$par, 
                    nrow = years.samp, ncol = 3, byrow = T)
dimnames(diff_pars) <- list(c(2015:(2014+years.samp)), c('ju', 'subad', 'ad'))

out_pars
diff_pars