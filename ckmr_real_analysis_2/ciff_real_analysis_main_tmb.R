# =================================================================
# CIFFs actual analysis of real data - ROUND 2 
# Author: Luke Lloyd-Jones and Isabel Peters
# Date started: 25/11/2021
# Date updated: 02/12/2021
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
compile("ciff_cpp/ciff_nllk.cpp")
compile("ciff_cpp/ciff_nllk_pmales.cpp")
dyn.load(dynlib("ciff_cpp/ciff_nllk"))
dyn.load(dynlib("ciff_cpp/ciff_nllk_pmales"))
lethal = 0
#------------------------------------------------------------------------
# Read back in the kin sets and the pop sample and prime the combinations
# and kin for fitting.
# All kin and combinations have been generated in previous setup in the
# ../kin
# ------------------------------------------------------------------------
mat.kin <- get(load("../kin_and_combos/comp_kin_out/mat_kin.Rdata"))
#  -----------------------------------------------------------------------
# Fit the CKMR non-lethal model
#  -----------------------------------------------------------------------
# -----------------
# Do the TMB setup
# -----------------
Yad_range <- as.numeric(colnames(mat.kin$mop.comp))
Bju_range <- as.numeric(colnames(mat.kin$mhsp.comp))

data <- list(comb_mop  = as.array(mat.kin$mop.comp),    
             kin_mop   = as.array(mat.kin$mop.kin),
             comb_mhsp = as.array(mat.kin$mhsp.comp),       
             kin_mhsp  = as.array(mat.kin$mhsp.kin),
             symin     = min(Yad_range), 
             symax     = max(Yad_range), 
             bymin     = min(Bju_range),  
             bymax     = max(Bju_range), 
             lethal    = 0, 
             fudge_factor =  1  - 0.081)
map  <- list(Z1 = factor(NA))
#map  <- list()
Z.val <- 0.6
ciff_nllk <- MakeADFun(data, 
                       parameters = list(N0_f = log(500), Z1 = Z.val, R1 = 0.0),
                       map = map,
                       DLL = "ciff_nllk")
ciff_nllk$control <- list(REPORT = 1)
# -----------
# Numeric fit
# -----------
fit1 = nlminb(ciff_nllk$par, 
              ciff_nllk$fn, 
              ciff_nllk$gr, 
              lower = c(1,    -20, -10), 
              upper = c(20000, 20, 10))
rep <- sdreport(ciff_nllk)
summary(rep)
print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk$gr(fit1$par))) 
print(fit1$par)
ciff_nllk$env$report()$N
fit1$par
# ---------------------------------------------------------------------------
# Given the fits, plot the abundance against the truth and write out a plot 
# to a results folder
# ---------------------------------------------------------------------------
years.samp <- 5
rep <- sdreport(ciff_nllk)
N.res <- as.data.frame(summary(rep)[grep("log.", rownames(summary(rep))), ])
# Build the results into a item to pass to the plot function
years <- seq(2013, 2019)
N.res$CV <- sqrt(exp(N.res[, 2]^2) - 1)
jpeg(paste0("plots/real_abundance_plot_Z", Z.val, ".jpeg"),
     pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 900, height = 800)
plot_abund_cv(years, N.res$Estimate, N.res$CV, N.res$"Std. Error", 0, 0, 7500)
dev.off()

# ===========================================
# Try again with females annnnnd males
# ===========================================

pat.kin <- get(load("../kin_and_combos/comp_kin_out/pat_kin.Rdata"))
Yad_range <- as.numeric(colnames(mat.kin$mop.comp))
Bju_range <- as.numeric(colnames(mat.kin$mhsp.comp))
diag(pat.kin$phsp.kin) <- 0
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
             surv_prior = 1,
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


# ---------------------------------------------------------------------------
# Summarise the tables generated
# ---------------------------------------------------------------------------
library(knitr)
# Abundance table
kable(round(abund.ests[-1, ], 0), format = "latex", booktabs = T)
kable(N.true, format = "latex", booktabs = T)
# Par ests table
thata.ests <- signif(abs(thata.ests), 3)
kable(thata.ests[-1, ], format = "latex", booktabs = T)
# ---------------------------------------------------------------------------
# Looking at the likelihood 
# 01/05/2020
# ---------------------------------------------------------------------------
com.kin$n_mhsps[2014, 2015] <- 2
com.kin$n_mhsps[2015, 2016] <- 0
data <- list(comb_mop  = as.array(com.kin$n_comp_MOP_noage), kin_mop = as.array(com.kin$n_MOP_noage),
             comb_mhsp = as.array(com.kin$n_comp_HSP), kin_mhsp = as.array(com.kin$n_mhsps),
             symin  = min(com.kin$Yad_range), symax = max(com.kin$Yad_range), 
             bymin = min(com.kin$Bju_range),  bymax = max(com.kin$Bju_range), 
             lethal = com.kin$lethal, fudge_factor =  1  - 0.1683315)
ciff_nllk <- MakeADFun(data, 
                       parameters = list(N0_f = log(500), Z1 = 0.3, R1 = 0.0),
                       map = map,
                       DLL = "ciff_nllk")
N0_f_vals <- seq(500, 12000, 10) 
llks      <- array(0, length(N0_f_vals))
for (i in seq(1, length(N0_f_vals)))
{
  llks[i] <- ciff_nllk$fn(log(N0_f_vals[i]))
}
plot(N0_f_vals, -llks)
llks_2 <- llks
llks_3 <- llks
llks_4 <- llks
llks_5 <- llks
llks_glob <- llks_2 - llks_3 - llks_4 - llks_5
plot(N0_f_vals, llks_glob, ylim = c(-20, 0))
lines(N0_f_vals, -llks_2)
lines(N0_f_vals, -llks_3)
lines(N0_f_vals, -llks_4)
lines(N0_f_vals, -llks_5)

round(ciff_nllk$report()$exp_mhsps)
# ---------------------------------
# Doing the Poisson thingo
# ---------------------------------
llks_delta <- list()
k <- 0
for (j in seq(0, 0.9, 0.1))
{
  k <- k + 1
  llks_pois      <- array(0, length(N0_f_vals))
  for (i in seq(1, length(N0_f_vals)))
  {
  ciff_nllk <- MakeADFun(data, 
                         parameters = list(N0_f = log(N0_f_vals[i]), Z1 = j, R1 = 0.0),
                         map = map,
                         DLL = "ciff_nllk")
  #exp_mhsps           <- ciff_nllk$report()$exp_mhsps
  exp_mhsps           <- ciff_nllk$report()$mhsp_prob * as.array(data$comb_mhsp)
  #exp_mops            <- ciff_nllk$report()$exp_mops
  exp_mops            <- ciff_nllk$report()$mop_prob * as.array(data$comb_mop)
  exp_matS_at_1       <- exp_mhsps[1, 2] + exp_mhsps[2, 3] + exp_mhsps[3, 4] + exp_mhsps[4, 5]
  exp_matS_at_2ormore <- sum(exp_mhsps[1, 3:5]) + sum(exp_mhsps[2, 4:5]) + exp_mhsps[3, 5]
  exp_MOPs            <- sum(exp_mops)
  llks_pois[i] <- sum( dpois(2:5, exp_matS_at_1, log = T) + dpois(0, exp_matS_at_2ormore, log = T) + dpois(0, exp_MOPs, log = T))
  }
  llks_delta[[k]]      <- llks_pois   
  names(llks_delta)[k] <- paste0(j)
  #plot(N0_f_vals, llks_pois)     
}     
plot(N0_f_vals, exp(llks_delta$`0`), type = "l", lwd = 3.5, 
     xlab = "No females", ylab = "likelihood", xlim = c(0, 12000),
     ylim = c(0, 1.4e-5))     
lines(N0_f_vals, exp(llks_delta$`0.1`), lwd = 3.5)#, col = 'green')
lines(N0_f_vals, exp(llks_delta$`0.2`), lwd = 3.5)#, col = 'red')     
lines(N0_f_vals, exp(llks_delta$`0.3`), lwd = 3.5)#, col = 'blue')  
lines(N0_f_vals, exp(llks_delta$`0.4`), lwd = 3.5)#, col = 'orange') 
lines(N0_f_vals, exp(llks_delta$`0.5`), lwd = 3.5)#, col = 'purple') 
lines(N0_f_vals, exp(llks_delta$`0.6`), lwd = 3.5)#, col = 'darkgreen') 
lines(N0_f_vals, exp(llks_delta$`0.7`), lwd = 3.5)#, col = 'yellow') 
lines(N0_f_vals, exp(llks_delta$`0.8`), lwd = 3.5)#, col = 'black') 
lines(N0_f_vals, exp(llks_delta$`0.9`), lwd = 3.5)#, col = 'cyan') 

plot(N0_f_vals, llks_delta$`0`, type = "l", lwd = 3.5, 
     xlab = "No females", ylab = "likelihood", xlim = c(0, 12000))     
lines(N0_f_vals, llks_delta$`0.1`, lwd = 3.5, col = 'green')
lines(N0_f_vals, llks_delta$`0.2`, lwd = 3.5, col = 'red')     
lines(N0_f_vals, llks_delta$`0.3`, lwd = 3.5, col = 'blue')  
lines(N0_f_vals, llks_delta$`0.4`, lwd = 3.5, col = 'orange') 
lines(N0_f_vals, llks_delta$`0.5`, lwd = 3.5, col = 'purple') 
lines(N0_f_vals, llks_delta$`0.6`, lwd = 3.5, col = 'darkgreen') 
lines(N0_f_vals, llks_delta$`0.7`, lwd = 3.5, col = 'yellow') 
lines(N0_f_vals, llks_delta$`0.8`, lwd = 3.5, col = 'black') 
lines(N0_f_vals, llks_delta$`0.9`, lwd = 3.5, col = 'cyan') 

# An example
mopsprobs <- round(mop.prob, 4)
mopsprobs[4, 2] <- 0
colnames(mopsprobs) <- colnames(as.array(data$comb_mop))
rownames(mopsprobs) <- rownames( as.array(data$comb_mop))
round(exp_mops, 2)

exp_mhsps     <- ciff_nllk$report()$mhsp_prob * as.array(data$comb_mhsp)
mhsps.prob.2  <- round(ciff_nllk$report()$mhsp_prob, 4)
mhsps.prob.2[2, 1] <- 0
diag(mhsps.prob.2) <- 0
colnames(mhsps.prob.2) <- colnames(as.array(data$comb_mhsp))
rownames(mhsps.prob.2) <- rownames( as.array(data$comb_mhsp))
round(exp_mhsps, 2)
