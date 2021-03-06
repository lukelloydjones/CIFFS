# =================================================================
# CIFFs simulation main
# Author: Luke Lloyd-Jones and Isabel Peters
# Date started: 04/10/2019
# Date updated: 06/12/2019
# =================================================================

# Libraries needed
suppressPackageStartupMessages({
  library(atease)   
  library(mvbutils)
  library(dplyr)
  library(offarray)
  library(reshape)
  library(RColorBrewer)
  library(TMB)
})

# Source functions required for simulation and analysis
source("functions/function_list.R")

# TMB compilation and linking
compile("ciff_cpp/ciff_nllk.cpp")
dyn.load(dynlib("ciff_cpp/ciff_nllk"))

# Some initial directories
# create possibly multiple directories if they do not currently exist.
make_directories <- function(...) {
  dirs <- unlist(list(...))
  stopifnot(length(dirs) > 0, is.character(dirs))
  for(d in dirs) {
    if(!dir.exists(d)) dir.create(d)
  }
  invisible(d)
}

# Some initial directories
make_directories("populations", 
                 "sampled_populations",
                 "found_kin", 
                 "preped_data", 
                 "results")

# ---------------------
# Simulate populations
#  - Creates three populations with different pop size trajectories 
#    in the 'population' subdirectory
#  - Warning: creates data in directories but only very small (kb)
# ---------------------

# Set the scenarios
mort.pars <- list(0.3) #, 0.23, 0.258, 0.30) # Increasing, stable, decreasing
# Funtion args population start size in 1978, generations and mortality parameter
start.size <- 10000
gens       <- 40
lapply(mort.pars, function(x) ciff_simulate(start.size, gens, x)) 

# ---------------------
# Sample the population
# ---------------------
sim.no        <- 30   # Number of simulations 
years.samp    <- 5    # How many years do you want to sample
sam.no.py     <- 100  # How many individuals do you want to sample per year
pop.samp.dead <- 0    # 1 => Lethal sampling, 0 => mark-recapture
scen          <- 0.3  # Of the scenarios simulated, which one do you want to analayse?
make_directories(file.path("sampled_populations", paste0("mort_", scen)))
rep.seq <- array(seq(1, sim.no))
# Sample the populations. Warning: write data to subdirectories
apply(rep.seq, 1, function(x) pop_sample(sam.no.py, years.samp, scen, gens, pop.samp.dead, x))

# ------------------------------------------------------------------------
# Read back in the samples and compute the kin. Save each as an Rdata file
# ------------------------------------------------------------------------
make_directories(file.path("found_kin", paste0("mort_", scen)))
apply(rep.seq, 1, function(x) find_kin(sam.no.py, years.samp, scen, x))

# ------------------------------------------------------------------------
# Read back in the kin sets and the pop sample and prime the combinations
# and kin for fitting
# ------------------------------------------------------------------------
make_directories(file.path("preped_data", paste0("mort_", scen)))
apply(rep.seq, 1, function(x) data_prep(sam.no.py, years.samp, scen, x, pop.samp.dead))

#  -----------------------------------------------------------------------
# Fit the CLMR lethal or non-lethal model
#  -----------------------------------------------------------------------

par.ests <- matrix(0, nrow = sim.no, ncol = 6) # Alwats a six parameter model
abundance.ests2 <- matrix(0, nrow = sim.no, ncol = years.samp + 2)
kin.nos <- matrix(0, nrow = length(rep.seq), ncol = 2)
for (repi in rep.seq)
{
  print(paste0("Fitting sample ", repi))
  load(paste0("preped_data/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
              years.samp, "_samp_no_py_", sam.no.py, ".Rdata"))
  # -----------------
  # Do the TMB setup
  # -----------------
  data  <- list(comb_mop = as.array(com.kin$n_comp_MOP_noage),  
                kin_mop = as.array(com.kin$n_MOP_noage),
                comb_fop = as.array(com.kin$n_comp_FOP_noage),  
                kin_fop = as.array(com.kin$n_FOP_noage),
                comb_mhsp = as.array(com.kin$n_comp_HSP), 
                kin_mhsp = as.array(com.kin$n_mhsps),
                comb_fhsp = as.array(com.kin$n_comp_HSP), 
                kin_fhsp = as.array(com.kin$n_phsps),
                symin  = min(com.kin$Yad_range), symax = max(com.kin$Yad_range), 
                bymin = min(com.kin$Bju_range),  bymax = max(com.kin$Bju_range), 
                lethal = com.kin$lethal)
  ciff_nllk <- MakeADFun(data, 
                         parameters = list(N0_f = log(1000), 
                                           N0_m = log(1000), 
                                           Z1 = 0.4, 
                                           Z2 = 0.45, R1 = 0.05, R2 = 0.01),
                         DLL = "ciff_nllk")
  ciff_nllk$control <- list(REPORT=1)
  # -----------
  # Numeric fit
  # -----------
  fit1 = nlminb(ciff_nllk$par, ciff_nllk$fn, ciff_nllk$gr, 
                lower = c(1, 1, -10, -12, -10, -10), 
                upper = c(20000, 20000, 10, 10, 10, 10))
  sdreport(ciff_nllk, fit1$par)
  print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk$gr(fit1$par))) 
  print(fit1$par)
  abundance.ests2[repi, ] <- ciff_nllk$env$report()$N
  par.ests[repi, ] <- fit1$par
  print("Abundance through time")
  print(com.kin$N)
  print("No. POPs")
  print(sum(com.kin$n_MOP_noage) + sum(com.kin$n_FOP_noage))
  kin.nos[repi, 1] <- sum(com.kin$n_MOP_noage) + sum(com.kin$n_FOP_noage)
  print("No. HSPs")
  kin.nos[repi, 2] <- sum(com.kin$n_mhsps) + sum(com.kin$n_phsps)
  print(sum(com.kin$n_mhsps) + sum(com.kin$n_phsps))
  #Sys.sleep(2)
}

# ---------------------------------------------------------------------------
# Given the fits, plot the abundance against the truth and write out a plot 
# to a results folder
# ---------------------------------------------------------------------------
colnames(abundance.ests2) <- com.kin$years
plot_abundance(abundance.ests2, sam.no.py, years.samp, gens, scen, pop.samp.dead)
colMeans(par.ests)
colMeans(abundance.ests2)
apply(abundance.ests2, 2, function(x) median(x))
# Standard errors
sds.N    <- apply(abundance.ests2, 2, function(x) sd(x))
sds.pars <- apply(par.ests, 2, function(x) sd(x))
# CVs
sds.N / colMeans(abundance.ests2)
sds.pars / colMeans(par.ests)
# Expected number of kin
colMeans(kin.nos)
