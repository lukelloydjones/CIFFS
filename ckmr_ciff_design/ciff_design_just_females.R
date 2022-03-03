# =============================================================================
# Design update to really try to get at the data. Males are no good now
# so we try and focus on the just the females. Can we do it? We assume we
# are just going to go for 728 individuals
# Author: Luke Lloyd-Jones
# Date: 07/04/2020
# Date updated: 26/05/2020
#   7/05/2020 added a part to look likelihood for different rho, delta, N_f
#   26/05/2020 changed mortality to 0.35, 0.5, 0.7  such that it's more like 
#              the survival range of interest
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
#compile("ciff_cpp/ciff_prob_deriv_fem.cpp")
dyn.load(dynlib("ciff_cpp/ciff_prob_deriv_fem"))

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

ckmr_design_func <- function(theta, type, mort, fem_size, roi, max_total_samps = 0, 
                             draw_plot = 0, min = 0, max = 5000)
{
  # print(prob_vec)
  # print(type)
  # PARS
  #  theta[1:years.samp] - vector of proportion of individuals to sample in each of 1:years.samp years
  #  theta[years.samp+1] - prop_ju_py - proportion of juveniles in each year
  #  theta[years.samp+2] - prop_sa_py - proportion of aub-adults in each year
  #  theta[years.samp+3] - prop_ad_py - proportion of adulds in each year
  #  type                - what type of objective function to return 1 -min CV, 2 - sum CV over years, 3 the CV vector for plotting
  #  mort                - mortality rate
  #  fem_size            - assumed initial female population size
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
    #pars <- c(log(800), log(874), 0.32, 0.29, -0.045, -0.050) # Hypothesised truth
    pars <- c(log(fem_size), mort, roi) # Hypothesised truth
    #800, 874, 0.32, 0.29, -0.045, -0.050
    #600, 600, 1.79, 1.79, -1,-1
    print("PARS")
    print(pars)
    N0_f <- pars[1]
    Z1   <- pars[2]
    R1   <- pars[3]
    
    environment(calc_ckmr_hess_fem) <- com.kin
    ckmr_hess <- calc_ckmr_hess_fem(pars)
    
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
                                   parameters = list(N0_f=N0_f, Z1=Z1, R1=R1), silent=TRUE)
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
    mhsps <- ciff_abund_grad$report()$mhsp_prob * com.kin$n_comp_HSP 
    mops <- ciff_abund_grad$report()$mop_prob * com.kin$n_comp_MOP_noage 
    
    # Return different objective functions
    print("****************************************")
    print("Samps py")
    print(samp.no.py)
    print("MIN CV, SUM OF CVs, CV vector")
    print(c(min(abund_cv), sum(abund_cv)))
    print(abund_cv)
    print("Expected number of MOPs and MHSPs")
    print(c(sum(mops, na.rm = T), sum(mhsps, na.rm = T)))
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
    } else if (type == 3) 
    {
      return(list(min(abund_cv), c(sum(mops, na.rm = T), sum(mhsps, na.rm = T))))
    } else {
      return(list(mhsps, mops, com.kin$n_comp_HSP, com.kin$n_comp_MOP_noage))
    }
    
    #return(sum(abund_cv) + diff)
    #return(list(min(abund_cv), sum(abund_cv), abund_cv)) # Minimum over CV vector, Sum of CV vector
    
  }, error=function(e) {
    print(paste0("ERROR - ", geterrmessage()))
    return(100)
  })
}

# ------------------------------------------------------
# Use all 728 inds and female proportion of 0.35
# ------------------------------------------------------
design_pars <- data.frame(rep(c(500, 1000, 2000), each = 12 ), 
                          rep(c(-0.1, -0.05, 0, 0.05), times = 9),
                          rep(rep(c(0.35, 0.5, 0.7), each = 4), times = 3))
colnames(design_pars) <- c("N_FEM", "ROI", "MORT")
min.cvs <- matrix(0, nrow = dim(design_pars)[1])
exp.kin <- matrix(0, nrow = dim(design_pars)[1], ncol = 3)
for (i in seq(1, dim(design_pars)[1]))
{
  des.res.i <- ckmr_design_func(c(68, 10, 27, 
                                  115, 32, 56, 
                                  77, 26, 55, 
                                  79, 2, 31, 
                                  77, 9, 64), 3, 
                                  fem_size = design_pars[i, 1], mort = design_pars[i, 3], roi = design_pars[i, 2], 
                                  draw_plot = 1, min = 0, max = 4000)
  min.cvs[i, ] <- des.res.i[[1]][1]
  exp.kin[i, 1] <- des.res.i[[2]][1]
  exp.kin[i, 2] <- des.res.i[[2]][2]
  exp.kin[i, 3] <- exp.kin[i, 1] + exp.kin[i, 2]
}
results <- cbind(design_pars, min.cvs, exp.kin)
results[, c("1", "2", "3")]  <- round(results[, c("1", "2", "3")], 0)
results[, 4] <- round(results[, 4], 3)*100
#results[, c(2, 3)] <- results[, c(3, 2)]
results$CV <- results$min.cvs
results <- results[, -4]
results$SURV = round(exp(-results$MORT), 2)
results <- results[, c(1, 3, 8, 2, 4, 5, 6, 7)]
# Write out for TEx
library(knitr)
kable(results, format = 'latex')

# ------------------------------------------------------
# Just look three scenarios 
# ------------------------------------------------------
design_pars <- data.frame(c(500, 1000, 1500, 2000), 
                          c(0, 0, 0, 0),
                          rep(0.5, each = 4))
colnames(design_pars) <- c("N_FEM", "ROI", "MORT")
min.cvs <- matrix(0, nrow = dim(design_pars)[1])
exp.kin <- matrix(0, nrow = dim(design_pars)[1], ncol = 3)
# Scenario 1
i=1
des.res.i <- ckmr_design_func(c(68, 10, 27, 
                                115, 32, 56, 
                                77, 26, 55, 
                                79, 2, 31, 
                                77, 9, 64), 4, 
                                fem_size = design_pars[i, 1], mort = design_pars[i, 3], roi = design_pars[i, 2], 
                                draw_plot = 1, min = 0, max = 4000)
save(des.res.i, file = "design_scenario_1.Rdata")
# -----
# MHSPs
# -----
hsp_kin  <- as.array(des.res.i[[1]])
hsp_comb <- as.array(des.res.i[[3]])
hsp_kin_chrts  <- ceiling(sapply(split(hsp_kin,  col(hsp_kin) - row(hsp_kin)), sum))
hsp_comb_chrts <- ceiling(sapply(split(hsp_comb, col(hsp_comb) - row(hsp_comb)), sum))
hsp_kin_chrts  <- hsp_kin_chrts[-seq(1, 6)]
hsp_comb_chrts <- hsp_comb_chrts[-seq(1, 6)]
# ----
# MOPs
# ----
mop_kin  <- as.array(des.res.i[[2]])
mop_comb <- as.array(des.res.i[[4]])
mop_kin_chrts  <- ceiling(sapply(split(mop_kin,  col(mop_kin)  - row(mop_kin)),  sum))
mop_comb_chrts <- ceiling(sapply(split(mop_comb, col(mop_comb) - row(mop_comb)), sum))

mop_comb_chrts[5] <- mop_comb_chrts[5] + mop_comb_chrts[3]
mop_comb_chrts[6] <- mop_comb_chrts[6] + mop_comb_chrts[2]
mop_comb_chrts[7] <- mop_comb_chrts[7] + mop_comb_chrts[1]
mop_comb_chrts    <- mop_comb_chrts[-seq(1, 4)]

mop_kin_chrts[5] <- mop_kin_chrts[5] + mop_kin_chrts[3]
mop_kin_chrts[6] <- mop_kin_chrts[6] + mop_kin_chrts[2]
mop_kin_chrts[7] <- mop_kin_chrts[7] + mop_kin_chrts[1]
mop_kin_chrts    <- mop_kin_chrts[-seq(1, 4)]
# ==========
# Scenario 2
# ==========
i=2
des.res.i <- ckmr_design_func(c(68, 10, 27, 
                                115, 32, 56, 
                                77, 26, 55, 
                                79, 2, 31, 
                                77, 9, 64), 4, 
                              fem_size = design_pars[i, 1], mort = design_pars[i, 3], roi = design_pars[i, 2], 
                              draw_plot = 1, min = 0, max = 4000)
save(des.res.i, file = "design_scenario_2.Rdata")
# -----
# MHSPs
# -----
hsp_kin  <- as.array(des.res.i[[1]])
hsp_comb <- as.array(des.res.i[[3]])
hsp_kin_chrts  <- ceiling(sapply(split(hsp_kin,  col(hsp_kin)  - row(hsp_kin)),  sum))
hsp_comb_chrts <- ceiling(sapply(split(hsp_comb, col(hsp_comb) - row(hsp_comb)), sum))
hsp_kin_chrts  <- hsp_kin_chrts[-seq(1, 6)]
hsp_comb_chrts <- hsp_comb_chrts[-seq(1, 6)]
# ----
# MOPs
# ----
mop_kin  <- as.array(des.res.i[[2]])
mop_comb <- as.array(des.res.i[[4]])
mop_kin_chrts  <- ceiling(sapply(split(mop_kin,  col(mop_kin)  - row(mop_kin)),  sum))
mop_comb_chrts <- ceiling(sapply(split(mop_comb, col(mop_comb) - row(mop_comb)), sum))

mop_comb_chrts[5] <- mop_comb_chrts[5] + mop_comb_chrts[3]
mop_comb_chrts[6] <- mop_comb_chrts[6] + mop_comb_chrts[2]
mop_comb_chrts[7] <- mop_comb_chrts[7] + mop_comb_chrts[1]
mop_comb_chrts    <- mop_comb_chrts[-seq(1, 4)]

mop_kin_chrts[5] <- mop_kin_chrts[5] + mop_kin_chrts[3]
mop_kin_chrts[6] <- mop_kin_chrts[6] + mop_kin_chrts[2]
mop_kin_chrts[7] <- mop_kin_chrts[7] + mop_kin_chrts[1]
mop_kin_chrts    <- mop_kin_chrts[-seq(1, 4)]
# ==========
# Scenario 3
# ==========
i=3
des.res.i <- ckmr_design_func(c(68, 10, 27, 
                                115, 32, 56, 
                                77, 26, 55, 
                                79, 2, 31, 
                                77, 9, 64), 4, 
                              fem_size = design_pars[i, 1], mort = design_pars[i, 3], roi = design_pars[i, 2], 
                              draw_plot = 1, min = 0, max = 4000)
save(des.res.i, file = "design_scenario_3.Rdata")
# -----
# MHSPs
# -----
hsp_kin  <- as.array(des.res.i[[1]])
hsp_comb <- as.array(des.res.i[[3]])
hsp_kin_chrts  <- ceiling(sapply(split(hsp_kin,  col(hsp_kin)  - row(hsp_kin)),  sum))
hsp_comb_chrts <- ceiling(sapply(split(hsp_comb, col(hsp_comb) - row(hsp_comb)), sum))
hsp_kin_chrts  <- hsp_kin_chrts[-seq(1, 6)]
hsp_comb_chrts <- hsp_comb_chrts[-seq(1, 6)]
# ----
# MOPs
# ----
mop_kin  <- as.array(des.res.i[[2]])
mop_comb <- as.array(des.res.i[[4]])
mop_kin_chrts  <- ceiling(sapply(split(mop_kin,  col(mop_kin)  - row(mop_kin)),  sum))
mop_comb_chrts <- ceiling(sapply(split(mop_comb, col(mop_comb) - row(mop_comb)), sum))

mop_comb_chrts[5] <- mop_comb_chrts[5] + mop_comb_chrts[3]
mop_comb_chrts[6] <- mop_comb_chrts[6] + mop_comb_chrts[2]
mop_comb_chrts[7] <- mop_comb_chrts[7] + mop_comb_chrts[1]
mop_comb_chrts    <- mop_comb_chrts[-seq(1, 4)]

mop_kin_chrts[5] <- mop_kin_chrts[5] + mop_kin_chrts[3]
mop_kin_chrts[6] <- mop_kin_chrts[6] + mop_kin_chrts[2]
mop_kin_chrts[7] <- mop_kin_chrts[7] + mop_kin_chrts[1]
mop_kin_chrts    <- mop_kin_chrts[-seq(1, 4)]

mop_kin_chrts
hsp_kin_chrts


# ==========
# Scenario 4
# ==========
i=4
des.res.i <- ckmr_design_func(c(68, 10, 27, 
                                115, 32, 56, 
                                77, 26, 55, 
                                79, 2, 31, 
                                77, 9, 64), 4, 
                              fem_size = design_pars[i, 1], mort = design_pars[i, 3], roi = design_pars[i, 2], 
                              draw_plot = 1, min = 0, max = 4000)
save(des.res.i, file = "design_scenario_4.Rdata")





# ------------------------------------------------------
# A test for the data at 240
# ------------------------------------------------------
design_pars <- data.frame(c(500, 1000, 10000), 
                          rep(0, times = 3),
                          rep(0.3, each = 3))
colnames(design_pars) <- c("N_FEM", "ROI", "MORT")
min.cvs <- matrix(0, nrow = dim(design_pars)[1])
exp.kin <- matrix(0, nrow = dim(design_pars)[1], ncol = 3)
# Scenario 1
years.samp <- 3    # How many years do you want to sample
des.res.i <- ckmr_design_func(c(33, 3, 8, 
                                85, 22, 38, 
                                11, 17, 12), 4, 
                              fem_size = design_pars[i, 1], mort = design_pars[i, 3], roi = design_pars[i, 2], 
                              draw_plot = 1, min = 0, max = 4000)