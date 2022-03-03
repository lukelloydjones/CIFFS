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
  #library(ggplot2)
  #library(reshape)
  #library(RColorBrewer)
  library(TMB)
  #library(DEoptim)
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
years.samp    <- 5    # How many years do you want to sample

# ---------------------------------
# Now make a function
# ---------------------------------
# Example theta = c(c(90, 90, 90, 90, 90, 50, 50, 50, 50, 50), c(0.3, 0.2, 0.5))
# theta, prob_vec,
# probs, theta,

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
  # Example samp.no.py     <- c(90, 90, 90, 90, 90, 50, 50, 50, 50, 50)
  # roofie    <- rep(130, years.samp) * theta[1:years.samp]
  # samp.no.py <- floor(roofie)
  # diff      <- sum(roofie - samp.no.py)
    
  #Start set for testing
  probs <- (1/1000) * theta[years.samp+1:2]
  if (sum(probs) > 1)
  {
    stop("Probs out of bounds")
  }
  prob_vec = c(1-sum(probs), probs)
    
  samp.no.py <- floor(theta[1:years.samp])
  if (max_total_samps != 0 && sum(samp.no.py) > max_total_samps)
  {
    stop("Max total samps exceeded")
  }
  
  #rep(c(90, 50), 5)  # How many individuals do you want to sample per year
  lethal        <- 0    # 1 => Lethal sampling, 0 => mark-recapture
  scen          <- 0.3  # Of the scenarios simulated, which one do you want to analayse?

  prop_ju_py    <- prob_vec[1]
  prop_sa_py    <- prob_vec[2]
  prop_ad_py    <- prob_vec[3]
  
  # A - 84,  J - 28,  SA - 18
  #com.kin <- data_prep_design(samp.no.py, years.samp, prop_ju_py, prop_sa_py, prop_ad_py, lethal)
  
  f_prop = 0.45
  com.kin <- prep_comb_kin_2(samp.no.py, years.samp, prop_ju_py, prop_sa_py, prop_ad_py, f_prop, lethal)
  
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
    data  <- list(symin   = min(com.kin$Yad_range), symax = max(com.kin$Yad_range),
                  bymin   = min(com.kin$Bju_range), bymax = max(com.kin$Bju_range),
                  lethal  = com.kin$lethal,
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
  print("THETA")
  print(theta[1:years.samp])
  #print(theta[years.samp+1:3])
  print(prob_vec)
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

ckmr_design_func(c(rep(100,5), 300, 300), 2)
ckmr_design_func(c(rep(50,5), 532, 30), 1)





nlminb(theta, ckmr_design_func, gradient = NULL, hessian = NULL, age.s1, 2,
       control = list(step.min = 100),
       lower = c(rep(0.01, 10), c(0.01, 0.01, 0.01)),
       upper = c(rep(1, 10),    c(1, 1, 1)))


d <- matrix(c(rep(c(1, 200), 5), rep(c(0, 1000), 2)), nrow = 7, ncol = 2, byrow = T)
out_genoud <- genoud(ckmr_design_func, nvars=7, pop.size=60, max.generations=80, wait.generations=15,
                     hard.generation.limit=F,
                     Domains=d,solution.tolerance=1e-4,
                     boundary.enforcement=1, 
                     data.type.int=T, 
                     type = 1, max_total_samps = 500)


# OK what to do now
# -----------------
# Want to change sample size, years sampled, number per year and age class to see what we get


# Real data 
new.data <- read.csv("C:/Users/PET247/Desktop/CIFF_capturedata_geneticsamples.csv")
head(new.data)
new.data.wb <- new.data[which(new.data$Wing.Biopsy == "Y" | new.data$Wing.Biopsy == "y"), ]

table(new.data.wb$Sex)

#  F   M 
#  110 162   

table(format(as.Date(new.data.wb$Date, "%d/%m/%Y"), "%Y"))

# 2018 2019 
# 114  158

table(new.data.wb$AgeClass)
new.data.wb$Year <- format(as.Date(new.data.wb$Date, "%d/%m/%Y"), "%Y")
table(new.data.wb$AgeClass, new.data.wb$Year)

dim(new.data.wb)/96
17 +  46

# Back of the envelope 
# 268 + 240 = 508 total sample for 96 well plates that's six plates 5.33. After quality going to be at  
# The optimiser says take as many per year and 0.58 juveniles to 42 sub-adults. No adults.
# Cost for DartSeq too much. 1 plate DaRT seq 5000. Enzyme prep 5000. CartCap per plate 1800
# What's the minimum CV for the data given different scenarios of cost under assumption that kin finding will work.


