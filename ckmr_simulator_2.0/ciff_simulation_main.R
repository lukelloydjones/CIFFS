# =================================================================
# CIFFs simulation main
# Author: Luke Lloyd-Jones
# Date started: 04/10/2019
# Date updated: 16/10/2019
# =================================================================
# install.packages("/Users/PET247/Dropbox/smp_summer_isabel_peters/week_2/resources/offarray_1.0.105.tar.gz", repos = NULL, type="source")
# Libraries needed
suppressPackageStartupMessages({
library(atease)   
library(mvbutils)
#library(vecless)
library(dplyr)
library(offarray)
library(reshape)
library(RColorBrewer)
})
# Source functions required for simulation and analysis
source("functions/function_list.R")
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
mort.pars <- list(0.30) #, 0.23, 0.258, 0.30) # Increasing, stable, decreasing
# Funtion args population start size in 1978, generations and mortality parameter
start.size <- 10000
gens       <- 40
lapply(mort.pars, function(x) ciff_simulate(start.size, gens, x)) 
# ---------------------
# Sample the population
# ---------------------
sim.no        <- 30   # Number of simulations
years.samp    <- 5   # How many years do you want to sample
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
# ------------------------------------------------------- -----------------
make_directories(file.path("preped_data", paste0("mort_", scen)))
apply(rep.seq, 1, function(x) data_prep(sam.no.py, years.samp, scen, x, pop.samp.dead))
#  -----------------------------------------------------------------------
# Fit the CLMR lethal or non-lethal model
#  -----------------------------------------------------------------------
par.ests <- matrix(0, nrow = sim.no, ncol = 6) # Alwats a six parameter model
abundance.ests <- matrix(0, nrow = sim.no, ncol = years.samp + 2)
kin.nos <- matrix(0, nrow = length(rep.seq), ncol = 2)
for (repi in rep.seq)
{
  print(paste0("Fitting sample ", repi))
  load(paste0("preped_data/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
            years.samp, "_samp_no_py_", sam.no.py, ".Rdata"))
  environment(ciff_ckmr_loglik) <- com.kin
  # -----------
  # Numeric fit
  # -----------
  fit1 = nlminb(c(1000, 1000, 0.4, 0.4, 1.1, 1.1), ciff_ckmr_loglik,
                lower = c(1, 1, -10, -10, -10 , -11),
                upper = c(10000, 10000, 10, 10, 200, 200))
  print(fit1$par)
  abundance.ests[repi, ] <- com.kin$N
  par.ests[repi, ] <- fit1$par
  print("Abundance through time")
  print(com.kin$N)
  print("No. POPs")
  print(sum(com.kin$n_MOP_noage) + sum(com.kin$n_FOP_noage))
  kin.nos[repi, 1] <- sum(com.kin$n_MOP_noage) + sum(com.kin$n_FOP_noage)
  print("No. HSPs")
  kin.nos[repi, 2] <- sum(com.kin$n_mhsps) + sum(com.kin$n_phsps)
  print(sum(com.kin$n_mhsps) + sum(com.kin$n_phsps))
}
# ---------------------------------------------------------------------------
# Given the fits, plot the abundance against the truth and write out a plot 
# to a results folder
# ---------------------------------------------------------------------------
colnames(abundance.ests) <- com.kin$years+1
plot_abundance(abundance.ests, sam.no.py, years.samp, gens, scen, pop.samp.dead)
colMeans(par.ests)
colMeans(abundance.ests)
apply(abundance.ests, 2, function(x) median(x))
