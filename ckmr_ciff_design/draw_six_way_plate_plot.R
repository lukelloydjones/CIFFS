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


ckmr_design_func_pplot <- function(theta, type, max_total_samps = 0, draw_plot = 0, min = 0, max = 5000)
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

    
    samp.no.py <- matrix(floor(theta), nrow = years.samp, ncol = 3, byrow = T)
    if (max_total_samps != 0 && sum(samp.no.py) > max_total_samps)
    {
      stop("Max total samps exceeded")
    }
    dimnames(samp.no.py) <- list(c(2015:(2014+years.samp)), c('ju', 'subad', 'ad'))
    
    #rep(c(90, 50), 5)  # How many individuals do you want to sample per year
    lethal        <- 0    # 1 => Lethal sampling, 0 => mark-recapture
    scen          <- 0.3  # Of the scenarios simulated, which one do you want to analayse?
    
    f_prop = 0.45
    com.kin <- prep_comb_kin_2_var_props(samp.no.py, years.samp, f_prop, lethal)
    
    # ----------------------------------------------------------
    # Calculate the CKMR approximate to the LLGK hessian
    # ----------------------------------------------------------
    # Need a set of starting parameters Nf_0, Nm_0, Z1, Z2, R1, R2
    pars <- c(log(800), log(874), 0.32, 0.29, -0.045, -0.050) # Hypothesised truth
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
      #plot_abund_cv(com.kin$years, N, abund_cv, lethal=0)
      #plot_abund_cv(com.kin$years, N, abund_cv, lethal=0, 100, 3000)
      #plot_abund_cv(com.kin$years, N, abund_cv, lethal=0, 0, 3100)
      rep <- sdreport(ciff_abund_grad)
      N.res <- as.data.frame(summary(rep)[grep("log.", rownames(summary(rep))), 1])
      # Build the results into a item to pass to the plot function
      return(plot_abund_cv(com.kin$years, N.res, abund_cv, std.err, lethal = 0, min, max))

}

###########################################################
all     <-   ckmr_design_func_pplot(c(68, 10, 27, 
                                 115, 32, 56, 
                                 77, 26, 55, 
                                 79, 2, 31, 
                                 77, 9, 64), 1, draw_plot = 1, min = 500, max = 3000)
seven.p <-   ckmr_design_func_pplot(c(68,  10, 0, 
                                      115, 32, 27, 
                                      77,  26, 55, 
                                      79,  2, 31, 
                                      77,  9, 64), 1, draw_plot = 1, min = 500, max = 3000)
six.p   <-  ckmr_design_func_pplot(c(68,  10, 0, 
                                     115, 32, 0, 
                                     77,  26, 20, 
                                     79,  2, 20, 
                                     77,  9, 41), 1, draw_plot = 1, min = 500, max = 3000)
# 5 plates, 480 samples
five.p   <-  ckmr_design_func_pplot(c(63, 10,  0, 
                   105, 32, 0, 
                   77, 26, 0, 
                   79, 2, 0, 
                   77, 9, 0), 1, draw_plot = 1, min = 500, max = 5000)
# 4 plates, 384 samples
four.p   <-  ckmr_design_func_pplot(c(43,  10, 0, 
                   69,  32, 0, 
                   57,  26, 0, 
                   59,  2,  0, 
                   77,  9,  0), 1, draw_plot = 1, min = 500, max = 5000)
# 3 plates, 288 samples
three.p   <-  ckmr_design_func_pplot(c(20,  10, 0, 
                   30,  32, 0, 
                   37,  26, 0, 
                   49,  2,  0, 
                   73,  9,  0), 1, draw_plot = 1, min = 500, max = 5000)
# ------------------------------------------------
# You gotta plot it out plot it out
# ------------------------------------------------
jpeg(filename = "plate_od_six_way_simulations.jpg", pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1580, height = 1080)
multiplot(all, five.p, seven.p, four.p, six.p, three.p, cols = 3)
grid.text('Years modelled', x=.53, y=.013, gp=gpar(fontsize=30, face = "bold"))
grid.text('Adult abundance (count)', x=.01, y=.55, gp=gpar(fontsize=30, face = "bold"), rot = 90)
dev.off() 
# =========================================
# Quick tabling
# -------------
all     <-   data.frame(t(matrix(c(68, 10, 27, 
                                      115, 32, 56, 
                                      77, 26, 55, 
                                      79, 2, 31, 
                                      77, 9, 64), nrow = 3, ncol = 5)))
all$Total <- rowSums(all)
all  <- rbind(all, colSums(all))
rownames(all) <- c("2015", "2016", "2017", "2018", "2019", "Total")

seven.p <-   data.frame(t(matrix(c(68,  10, 0, 
                                      115, 32, 27, 
                                      77,  26, 55, 
                                      79,  2, 31, 
                                      77,  9, 64), nrow = 3, ncol = 5)))
seven.p$Total <- rowSums(seven.p)
seven.p  <- rbind(seven.p, colSums(seven.p))
rownames(seven.p) <- c("2015", "2016", "2017", "2018", "2019", "Total")

six.p   <-  data.frame(t(matrix(c(68,  10, 0, 
                                     115, 32, 0, 
                                     77,  26, 20, 
                                     79,  2, 20, 
                                     77,  9, 41), nrow = 3, ncol = 5)))
six.p$Total <- rowSums(six.p)
six.p  <- rbind(six.p, colSums(six.p))
rownames(six.p) <- c("2015", "2016", "2017", "2018", "2019", "Total")

# 5 plates, 480 samples
five.p   <- data.frame(t(matrix( c(63, 10,  0, 
                                      105, 32, 0, 
                                      77, 26, 0, 
                                      79, 2, 0, 
                                      77, 9, 0), nrow = 3, ncol = 5)))
five.p$Total <- rowSums(five.p)
five.p  <- rbind(five.p, colSums(five.p))
rownames(five.p) <- c("2015", "2016", "2017", "2018", "2019", "Total")

# 4 plates, 384 samples
four.p   <-  data.frame(t(matrix(c(43,  10, 0, 
                                      69,  32, 0, 
                                      57,  26, 0, 
                                      59,  2,  0, 
                                      77,  9,  0), nrow = 3, ncol = 5)))
four.p$Total <- rowSums(four.p)
four.p  <- rbind(four.p, colSums(four.p))
rownames(four.p) <- c("2015", "2016", "2017", "2018", "2019", "Total")

# 3 plates, 288 samples
three.p   <-  data.frame(t(matrix(c(20,  10, 0, 
                                       30,  32, 0, 
                                       37,  26, 0, 
                                       49,  2,  0, 
                                       73,  9,  0), nrow = 3, ncol = 5)))
three.p$Total <- rowSums(three.p)
three.p  <- rbind(three.p, colSums(three.p))
rownames(three.p) <- c("2015", "2016", "2017", "2018", "2019", "Total")


library(knitr)
kable(all, format = 'latex')
kable(seven.p, format = 'latex')
kable(six.p, format = 'latex')
kable(five.p, format = 'latex')
kable(four.p, format = 'latex')
kable(three.p, format = 'latex')
