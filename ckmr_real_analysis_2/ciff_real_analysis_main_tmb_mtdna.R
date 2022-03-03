# =================================================================
# CIFFs actual analysis of real data - ROUND 2 - WITH MTDNA ACTION
# Author: Luke Lloyd-Jones 
# Date started: 04/12/2021
# Date updated: 14/12/2021
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
  library(ggplot2)
})

# Source functions required for simulation and analysis
source("functions/function_list.R")
# TMB compilation and linking
compile("ciff_cpp/ciff_nllk_pmales.cpp")
compile("ciff_cpp/ciff_nllk_pmales_mtdna.cpp")
dyn.load(dynlib("ciff_cpp/ciff_nllk_pmales"))
dyn.load(dynlib("ciff_cpp/ciff_nllk_pmales_mtdna"))
#------------------------------------------------------------------------
# Read back in the kin sets and the pop sample and prime the combinations
# and kin for fitting.
# All kin and combinations have been generated in previous setup in the
# ../kin
# Loop over all possible comparisons for HSPs and just ask are they HSP
# but adding the PHSP and MHSP probabilities together.
# Then we have an additional loop over over those pairs for which the 
# additional data are available i.e., not for the missing ones which looks at 
# the conditional lglk of the pairwise datum given that the pair is kin of 
# whatever specific type (and the pair's other covariates).
# ------------------------------------------------------------------------
mat.kin <- get(load("../kin_and_combos/comp_kin_out/mat_kin.Rdata"))
pat.kin <- get(load("../kin_and_combos/comp_kin_out/pat_kin.Rdata"))

Yad_range <- as.numeric(colnames(mat.kin$mop.comp))
Bju_range <- as.numeric(colnames(mat.kin$mhsp.comp))
diag(pat.kin$phsp.kin)  <- 0
diag(pat.kin$phsp.comp) <- 0
data <- list(comb_mop  = as.array(mat.kin$mop.comp),    
             kin_mop   = as.array(mat.kin$mop.kin),
             comb_mhsp = as.array(mat.kin$mhsp.comp),       
             kin_mhsp  = as.array(mat.kin$mhsp.kin),
             comb_fop  = as.array(pat.kin$ppop.comp),    
             kin_fop   = as.array(pat.kin$ppop.kin),
             comb_fhsp = as.array(pat.kin$phsp.comp),       
             kin_fhsp  = as.array(pat.kin$phsp.kin),
             symin     = min(Yad_range), 
             symax     = max(Yad_range), 
             bymin     = min(Bju_range),  
             bymax     = max(Bju_range), 
             lethal    = 0, 
             surv_prior   = 1,
             fudge_factor =  1  - 0.081,
             Mu_Z1 = 0.62, 
             SD_Z1 = 0.2,
             Mu_Z2 = 0.48,
             SD_Z2 = 0.2)
map  <- list(Z1 = factor(NA))
map  <- list(Z1 = factor(NA), Z2 = factor(NA))
map  <- list()
Z1.val <- 0.6
Z2.val <- 0.48
ciff_nllk_pmal <- MakeADFun(data, 
                            parameters = list(N0_f = log(500), N0_m = log(1000), 
                                              Z1 = Z1.val, Z2 = Z2.val, 
                                              R1 = 0.0, R2 = 0.0),
                            map = map,
                            DLL = "ciff_nllk_pmales")
ciff_nllk_pmal$control <- list(REPORT = 1)
# -----------
# Numeric fit
# -----------

fit2 = nlminb(ciff_nllk_pmal$par, ciff_nllk_pmal$fn, ciff_nllk_pmal$gr, 
              lower = c(1,     1,     0.01, 0.01, -10, -10), 
              upper = c(20000, 20000, 1,     1,  10,  10))
rep <- sdreport(ciff_nllk_pmal)
summary(rep)
print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk_pmal$gr(fit2$par))) 
print(fit2$par)
ciff_nllk_pmal$env$report()$N_f

hist(rnorm(10000, 0.5, 0.15))

# ---------------------------------------------------------------------------
# Given the fits, plot the abundance against the truth and write out a plot 
# to a results folder
# ---------------------------------------------------------------------------
years.samp <- 5
rep <- sdreport(ciff_nllk_pmal)
N.res <- as.data.frame(summary(rep)[grep("log\\(N_f\\)", rownames(summary(rep))), ])
# Build the results into a item to pass to the plot function
years <- seq(2013, 2019)
N.res$CV <- sqrt(exp(N.res[, 2]^2) - 1)
jpeg(paste0("plots/real_female_abundance_plot_priorZ", Z.val, ".jpeg"),
     pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1000, height = 800)
plot_abund_cv(years, N.res$Estimate, N.res$CV, N.res$"Std. Error", 0, 0, 7500)
dev.off()





# ---------------------------------------------------------------

# Take out the male free parameters and male them same as female

# ---------------------------------------------------------------

data <- list(comb_mop  = as.array(mat.kin$mop.comp),    
             kin_mop   = as.array(mat.kin$mop.kin),
             comb_mhsp = as.array(mat.kin$mhsp.comp),       
             kin_mhsp  = as.array(mat.kin$mhsp.kin),
             comb_fop  = as.array(pat.kin$ppop.comp),    
             kin_fop   = as.array(pat.kin$ppop.kin),
             comb_fhsp = as.array(pat.kin$phsp.comp),       
             kin_fhsp  = as.array(pat.kin$phsp.kin),
             symin     = min(Yad_range), 
             symax     = max(Yad_range), 
             bymin     = min(Bju_range),  
             bymax     = max(Bju_range), 
             lethal    = 0, 
             surv_prior = 1,
             fudge_factor =  1  - 0.081,
             Mu_Z1 = 0.62, 
             SD_Z1 = 0.2)

ciff_nllk_pmal <- MakeADFun(data, 
                            parameters = list(N0_f = log(500), 
                                              N0_m = log(1000), 
                                              Z1 = Z1.val,
                                              R1 = 0.0),
                            map = map,
                            DLL = "ciff_nllk_pmales_mtdna")
ciff_nllk_pmal$control <- list(REPORT = 1)
# -----------
# Numeric fit
# -----------

fit2 = nlminb(ciff_nllk_pmal$par, ciff_nllk_pmal$fn, ciff_nllk_pmal$gr, 
              lower = c(1,     1,     0.01, -10), 
              upper = c(20000, 20000, 1,     10))
rep <- sdreport(ciff_nllk_pmal)
summary(rep)
print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk_pmal$gr(fit2$par))) 


years.samp <- 5
rep <- sdreport(ciff_nllk_pmal)
N.res <- as.data.frame(summary(rep)[grep("log\\(N_f\\)", rownames(summary(rep))), ])
# Build the results into a item to pass to the plot function
years <- seq(2013, 2019)
N.res$CV <- sqrt(exp(N.res[, 2]^2) - 1)
jpeg(paste0("plots/real_female_abundance_plot_priorZ", Z.val, ".jpeg"),
     pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1000, height = 800)
plot_abund_cv(years, N.res$Estimate, N.res$CV, N.res$"Std. Error", 0, 0, 7500)
dev.off()

# ---------------------------------------------------------------

# Put the male as Nm_0 = Nf_0 * pi where 0 < pi <1

# ---------------------------------------------------------------


data <- list(comb_mop  = as.array(mat.kin$mop.comp),    
             kin_mop   = as.array(mat.kin$mop.kin),
             comb_mhsp = as.array(mat.kin$mhsp.comp),       
             kin_mhsp  = as.array(mat.kin$mhsp.kin),
             comb_fop  = as.array(pat.kin$ppop.comp),    
             kin_fop   = as.array(pat.kin$ppop.kin),
             comb_fhsp = as.array(pat.kin$phsp.comp),       
             kin_fhsp  = as.array(pat.kin$phsp.kin),
             symin     = min(Yad_range), 
             symax     = max(Yad_range), 
             bymin     = min(Bju_range),  
             bymax     = max(Bju_range), 
             lethal    = 0, 
             surv_prior = 1,
             fudge_factor =  1  - 0.081,
             Mu_Z1 = 0.62, 
             SD_Z1 = 0.2)

ciff_nllk_pmal <- MakeADFun(data, 
                            parameters = list(N0_f = log(500), 
                                              IP = log(0.8 / (1-0.8)), 
                                              Z1 = Z1.val,
                                              R1 = 0.0),
                            map = map,
                            DLL = "ciff_nllk_pmales_mtdna")
ciff_nllk_pmal$control <- list(REPORT = 1)
# -----------
# Numeric fit
# -----------

fit2 = nlminb(ciff_nllk_pmal$par, ciff_nllk_pmal$fn, ciff_nllk_pmal$gr, 
              lower = c(1,     -20,     0.01, -10), 
              upper = c(20000, 20, 1,     10))
rep <- sdreport(ciff_nllk_pmal)
summary(rep)
print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk_pmal$gr(fit2$par))) 




# ---------------------------------------------------------------

# Do the iteration now with coarse level loop
# Back here after mtDNA interlude 13/12/2021

# ---------------------------------------------------------------

compile("ciff_cpp/ciff_nllk_pmales_mtdna.cpp")
dyn.load(dynlib("ciff_cpp/ciff_nllk_pmales_mtdna"))

mop.kin <- get(load("../kin_and_combos/comp_kin_out/mop_kin.Rdata"))
fop.kin <- get(load("../kin_and_combos/comp_kin_out/fop_kin.Rdata"))
sib.kin <- get(load("../kin_and_combos/comp_kin_out/sibling_kin.Rdata"))

Yad_range <- as.numeric(colnames(mop.kin$mop.comp))
Bju_range <- as.numeric(colnames(sib.kin$sib.comp))

data <- list(comb_mop  = as.array(mop.kin$mop.comp),    
             kin_mop   = as.array(mop.kin$mop.kin),
             comb_fop  = as.array(fop.kin$ppop.comp),    
             kin_fop   = as.array(fop.kin$ppop.kin),
             comb_sib  = as.array(sib.kin$sib.comp),       
             kin_sib   = as.array(sib.kin$sib.kin),
             symin     = min(Yad_range), 
             symax     = max(Yad_range), 
             bymin     = min(Bju_range),  
             bymax     = max(Bju_range), 
             lethal    = 0, 
             surv_prior = 0,
             fudge_factor =  1  - 0.081,
             Mu_Z1 = 0.62, 
             SD_Z1 = 0.2)

Z1.val <- 0.3 # Just fix this for now so I can get the llgk working in R
map  <- list(Z1 = factor(NA))
#map  <- list()
ciff_nllk_pmal <- MakeADFun(data, 
                            parameters = list(N0_f = log(500), 
                                              IP = log(0.8 / (1 - 0.8)), 
                                              Z1 = Z1.val,
                                              R1 = 0.0),
                            map = map,
                            DLL = "ciff_nllk_pmales_mtdna")
ciff_nllk_pmal$control <- list(REPORT = 1)
ciff_nllk_pmal$fn()

# -----------
# Numeric fit
# -----------

fit2 = nlminb(ciff_nllk_pmal$par, ciff_nllk_pmal$fn, ciff_nllk_pmal$gr, 
              lower = c(1,     -20,     0.01, -10), 
              upper = c(20000, 20, 1,     10))
rep <- sdreport(ciff_nllk_pmal)
summary(rep)
print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk_pmal$gr(fit2$par))) 


years.samp <- 5
rep <- sdreport(ciff_nllk_pmal)
N.res <- as.data.frame(summary(rep)[grep("log\\(N_f\\)", rownames(summary(rep))), ])
# Build the results into a item to pass to the plot function
years <- seq(2013, 2019)
N.res$CV <- sqrt(exp(N.res[, 2]^2) - 1)
jpeg(paste0("plots/real_female_abundance_plot_priorZ", Z1.val, ".jpeg"),
     pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1000, height = 800)
plot_abund_cv(years, N.res$Estimate, N.res$CV, N.res$"Std. Error", 0, 0, 7500)
dev.off()

