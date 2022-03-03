# =================================================================
# CIFFs actual analysis of real data - ROUND 3 - WITH 
# even more MTDNA ACTION
# Author: Luke Lloyd-Jones 
# Date started: 16/12/2021
# Date updated: 17/01/2022
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

# ---------------------------------------------------------------

# Do the iteration now with coarse level loop
# Back here after mtDNA interlude 13/12/2021

# ---------------------------------------------------------------
# TMB compiling and linking

compile("ciff_cpp/ciff_nllk_pmales_mtdna.cpp")
dyn.load(dynlib("ciff_cpp/ciff_nllk_pmales_mtdna"))

# Read pre-computed kin matrices
mop.kin <- get(load("../kin_and_combos/comp_kin_out/mop_kin.Rdata"))
fop.kin <- get(load("../kin_and_combos/comp_kin_out/fop_kin.Rdata"))
sib.kin <- get(load("../kin_and_combos/comp_kin_out/sibling_kin.Rdata"))

# -----------------------------
# Read and setup the deffo sibs
# data matrix for TMB
# -----------------------------
defo.sibs <- read.csv("../mtDNA/results/all_siblings_plod34_mtdna_v_haplo.csv")
# Add a HSP and FSP index onto this as we need this in the clean up loop.
defo.sibs$FSP <- 0
defo.sibs$FSP[defo.sibs$PLOD > 100] <- 1
# Need to exclude those that have no haplo info
defo.sibs <- defo.sibs[-which(is.na(defo.sibs$HAP1) | 
                              is.na(defo.sibs$HAP2)), ]

# Need bc1, bc2, hap_frq2, same_hap, fsp_hsp
# The birth cohorts need to be adjust to represent the
# index of the matrix, which starts from 0 in C++
df.sibs.tmb <- defo.sibs[, c("IND1_BC", "IND2_BC", "HAP2_FRQ", "SAME_HAP", "FSP")]
df.sibs.tmb$IND1_BC <- df.sibs.tmb$IND1_BC - 2013
df.sibs.tmb$IND2_BC <- df.sibs.tmb$IND2_BC - 2013

Yad_range <- as.numeric(colnames(mop.kin$mop.comp))
Bju_range <- as.numeric(colnames(sib.kin$sib.comp))

# Some simple prior investigations
# X <- rnorm(100000, 0.57, 0.04)
# Y <- -log(X)

# Remove diagonal information

diag(sib.kin$sib.comp) <- 0
diag(sib.kin$sib.kin)  <- 0
df.sibs.tmb <- df.sibs.tmb[-which(df.sibs.tmb$IND1_BC == df.sibs.tmb$IND2_BC), ]

data <- list(comb_mop   = as.array(mop.kin$mop.comp),    
             kin_mop    = as.array(mop.kin$mop.kin),
             comb_fop   = as.array(fop.kin$ppop.comp),    
             kin_fop    = as.array(fop.kin$ppop.kin),
             comb_sib   = as.array(sib.kin$sib.comp),       
             kin_sib    = as.array(sib.kin$sib.kin),
             dsm_bc1    = as.array(df.sibs.tmb$IND1_BC),
             dsm_bc2    = as.array(df.sibs.tmb$IND2_BC),                      
             dsm_pi_h2  = as.array(df.sibs.tmb$HAP2_FRQ),                                      
             dsm_sm_hp  = as.array(df.sibs.tmb$SAME_HAP),                                             
             df_fsp_yn  = as.array(df.sibs.tmb$FSP),
             symin      = min(Yad_range), 
             symax      = max(Yad_range), 
             bymin      = min(Bju_range),  
             bymax      = max(Bju_range), 
             lethal     = 0, 
             surv_prior = 1,
             roi_prior  = 0,
             fng_rate   = 1  - 0.081,
             Mu_Z1      = 0.565, # 2 * log(0.57) - (1/2)*log(0.04^2 + 0.57^2), 2 * log(0.64) - (1/2)*log(0.02^2 + 0.64^2)  # 0.565
             SD_Z1      = 0.07,  # sqrt(-2 * log(0.57) + log(0.04^2 + 0.57^2)), sqrt(-2 * log(0.64) + log(0.15^2 + 0.64^2)) # 0.07
             Mu_R1      = 0,
             SD_R1      = 0.1)

map  <- list()
#map  <- list(R1 = factor(NA))
#map  <- list(Alpha = factor(NA))
ciff_nllk_all_sibs <- MakeADFun(data, 
                                parameters = list(N0_f  = log(500), 
                                                  IP    = log(0.8 / (1 - 0.8)), 
                                                  Z1    = 0.5,
                                                  R1    = 0,
                                                  Alpha = 0),
                                map = map,
                                DLL = "ciff_nllk_pmales_mtdna")
ciff_nllk_all_sibs$control <- list(REPORT = 1)
ciff_nllk_all_sibs$fn()

# -----------
# Numeric fit
# -----------

fit2 = nlminb(ciff_nllk_all_sibs$par, 
              ciff_nllk_all_sibs$fn, 
              ciff_nllk_all_sibs$gr, 
              lower = c(1,     -20, 0.01,  -10, -10, -10), 
              upper = c(20000,  20, 1.00,   10,  10,  20))
rep <- sdreport(ciff_nllk_all_sibs)
summary(rep)
print(paste("SCORE SCORE SCORE SCOOOOOOORRRREEEE", ciff_nllk_all_sibs$gr(fit2$par))) 
round(summary(rep), 3)

years.samp <- 5
N.res <- as.data.frame(summary(rep)[grep("log\\(N_f\\)", rownames(summary(rep))), ])
# Build the results into a item to pass to the plot function
years <- seq(2013, 2019)
N.res$CV <- sqrt(exp(N.res[, 2]^2) - 1)
#N.res$CV <- N.res$`Std. Error` / N.res$Estimate

jpeg(paste0("plots/real_female_abundance_plot_priorZ.jpeg"),
     pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1132.6, height = 700)
plot_abund_cv(years, N.res$Estimate, N.res$CV, N.res$"Std. Error", 0, 0, 3500, 0)
dev.off()


# ------------------------
# Some extras
# ------------------------
exp(7.020)*0.538

my <- exp(0.563) + 1/2*0.07^2
vy <- exp(2*0.563 + 0.07^2) * (exp(0.07^2) - 1)
sqrt(vy)
my
mean(y)
sd(y)
sd(y)^2

X <- rnorm(10000, 0.563, 0.07)
hist(X)

Y <- exp(-X)
hist(Y)
mean(Y)
sd(Y)
