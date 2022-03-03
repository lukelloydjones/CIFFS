# =============================================================================
# CODE to just do the final design plot
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
compile("ciff_cpp/ciff_prob_deriv_fem.cpp")
dyn.load(dynlib("ciff_cpp/ciff_prob_deriv_fem"))

source("functions/function_list.R")

# ---------------------------------
# Partition the data in myriad ways
# ---------------------------------
make_directories(file.path("preped_data"))
#ex.dir <- "/Users/llo080/Dropbox/CSIRO/summer_scholarship/smp_summer_isabel_peters/ckmr_simulator_2.0"
#ex.dir <- "/Users/pet247/Dropbox/smp_summer_isabel_peters/ckmr_simulator_2.0"
years.samp <- 5    # How many years do you want to sample



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
                                   parameters = list(N0_f=N0_f, Z1=Z1, R1=R1), silent=TRUE, DLL = "ciff_prob_deriv_fem")
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
    mops  <- ciff_abund_grad$report()$mop_prob * com.kin$n_comp_MOP_noage 
    
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
      return(plot_abund_cv(com.kin$years, N.res, abund_cv, std.err, lethal = 0, min, max))
    }
    
    # What type of objective function
    # if (type == 1)
    # {
    #   return(min(abund_cv))
    # } else if (type == 2)
    # {
    #   return(sum(abund_cv))
    # } else {
    #   return(list(min(abund_cv), c(sum(mops, na.rm = T), sum(mhsps, na.rm = T))))
    # }
    #return(plot_abund_cv(com.kin$years, N.res, abund_cv, std.err, lethal = 0, min, max))
    
    #return(sum(abund_cv) + diff)
    #return(list(min(abund_cv), sum(abund_cv), abund_cv)) # Minimum over CV vector, Sum of CV vector
    
  }, error=function(e) {
    print(paste0("ERROR - ", geterrmessage()))
    return(100)
  })
}

###########################################################
p1  <- ckmr_design_func(c(68, 10, 27, 
                                               115, 32, 56, 
                                               77, 26, 55, 
                                               79, 2, 31, 
                                               77, 9, 64), 3, 
                                             fem_size = 500, mort = 0.5, roi = -0.05, 
                                             draw_plot = 1, min = 0, max = 2500)
p2  <- ckmr_design_func(c(68, 10, 27, 
                           115, 32, 56, 
                           77, 26, 55, 
                           79, 2, 31, 
                           77, 9, 64), 3, 
                         fem_size = 1000, mort = 0.5, roi = -0.05, 
                         draw_plot = 1, min = 250, max = 2500)
p3  <- ckmr_design_func(c(68, 10, 27, 
                           115, 32, 56, 
                           77, 26, 55, 
                           79, 2, 31, 
                           77, 9, 64), 3, 
                         fem_size = 1000, mort = 0.5, roi = -0.1, 
                         draw_plot = 1, min = 250, max = 2500)
p4  <- ckmr_design_func(c(68, 10, 27, 
                           115, 32, 56, 
                           77, 26, 55, 
                           79, 2, 31, 
                           77, 9, 64), 3, 
                         fem_size = 2000, mort = 0.5, roi = -0.05, 
                         draw_plot = 1, min = 0, max = 6500)
# ------------------------------------------------
# You gotta plot it ou7t plot it out
# ------------------------------------------------
source("~/Dropbox/Post_Doc_QBI/R_Functions/ggplot_multi_func.R")
jpeg(filename = "design_four_way.jpg", pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1080, height = 780)
multiplot(p4, p3, p2, p1, cols = 2)
grid.text('Years modelled', x=.53, y=.018, gp=gpar(fontsize=30, face = "bold"))
grid.text('Female adult abundance (count)', x=.015, y=.56, gp=gpar(fontsize=30, face = "bold"), rot = 90)
dev.off() 

