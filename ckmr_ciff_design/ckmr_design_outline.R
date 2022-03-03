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
  library(reshape)
  library(RColorBrewer)
  library(TMB)
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
years.samp    <- 10    # How many years do you want to sample
sam.no.py     <- c(90, 90, 90, 90, 90, 50, 50, 50, 50, 50)
  #rep(c(90, 50), 5)  # How many individuals do you want to sample per year
lethal        <- 0    # 1 => Lethal sampling, 0 => mark-recapture
scen          <- 0.3  # Of the scenarios simulated, which one do you want to analayse?
repi          <- 1

# Start set for testing
prop_ju_py    <- 0.3
prop_sa_py    <- 0.2
prop_ad_py    <- 0.5
# A - 84,  J - 28,  SA - 18
com.kin <- data_prep_design(sam.no.py, years.samp, prop_ju_py, prop_sa_py, prop_ad_py, lethal)

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
ckmr_hess <- calc_ckmr_hess(pars, prop_ju_py, prop_sa_py, prop_ad_py)

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
                               parameters = list(N0_f=N0_f, N0_m=N0_m, Z1=Z1, Z2=Z2, R1=R1, R2=R2))
  abund_grad_year <- ciff_abund_grad$gr(pars)
  N <- ciff_abund_grad$report()$N
  names(N) <- as.character(com.kin$years)
  var             <- abund_grad_year %*% solve(ckmr_hess) %*% t(abund_grad_year)
  cv              <- sqrt(var) / N[as.character(year)]
  abund_cv[as.character(year)] = cv
}

# ---------------- 
# Draw a nice plot
# ----------------
plot_abund_cv(com.kin$years, N, abund_cv, lethal, min = 400, max = 2600)

# OK what to do now
# -----------------
# Want to change sample size, years sampled, number per year and age class to see what we get








