# =================================================================
# Small script to look at same versus different roost, kin and
# combinations for 
# Author: Luke Lloyd-Jones
# Date started: 20/01/2019
# Date updated: 02/04/2020
# =================================================================
#------------------------------------------------------------------------
# Read back in the kin sets and the pop sample and prime the combinations
# and kin for fitting
# ------------------------------------------------------------------------
kin.pops  <- read.table("for_analysis/pops_ids.txt")
kin.hsps  <- read.table("for_analysis/hsps_ids.txt")
#kin.hsps  <- read.table("for_analysis/hsps_pa_ids.txt")
ciff.samp <- read.table("for_analysis/ciff_set.txt", header = T)
big.set   <- read.csv("../mark_recapture_data/Copy of Copy of CIFF_capturedata_geneticsamples_tbu_CT .csv", header = T)
colnames(ciff.samp) <- c("id" ,        "Date"    ,   "mf"     , "Age",   "life.sg.status"   , "capture_loc", "year_samp", "Month"    , "year", "Number")
ciff.samp$mf[which(ciff.samp$mf == "F")] <- 0
ciff.samp$mf[which(ciff.samp$mf == "M")] <- 1
ciff.samp$life.sg.status[which(ciff.samp$life.sg.status == "Juvenile")] <- "J"
ciff.samp$life.sg.status[which(ciff.samp$life.sg.status == "Sub-adult")] <- "SA"
ciff.samp$life.sg.status[which(ciff.samp$life.sg.status == "Adult")] <- "A"
kin        <- list(kin.pops, kin.hsps)
names(kin) <- c("pops", "hsps")
com.kin    <- prep_comb_kin(kin, ciff.samp, lethal)
# -----------------------------
# First melt the ciff.samp data
# -----------------------------
ciff.samp.mlt <- data.frame(id_1    = rep(ciff.samp$id,       each = dim(ciff.samp)[1]), 
                            Date_1  = rep(ciff.samp$Date,     each = dim(ciff.samp)[1]), 
                            mf_1    = rep(ciff.samp$mf,       each = dim(ciff.samp)[1]), 
                            Age_1   = rep(ciff.samp$Age,      each = dim(ciff.samp)[1]), 
                            life.sg.status_1 = rep(ciff.samp$life.sg.status, each = dim(ciff.samp)[1]), 
                            capture_loc_1    = rep(ciff.samp$capture_loc,    each = dim(ciff.samp)[1]), 
                            year_samp_1      = rep(ciff.samp$year_samp,      each = dim(ciff.samp)[1]), 
                            Month_1          = rep(ciff.samp$Month,          each = dim(ciff.samp)[1]), 
                            year_1           = rep(ciff.samp$year,           each = dim(ciff.samp)[1]), 
                            Number_1         = rep(ciff.samp$Number,         each = dim(ciff.samp)[1]),
                            id_2    = rep(ciff.samp$id,             times = dim(ciff.samp)[1]), 
                            Date_2  = rep(ciff.samp$Date,           times = dim(ciff.samp)[1]), 
                            mf_2    = rep(ciff.samp$mf,             times = dim(ciff.samp)[1]), 
                            Age_2   = rep(ciff.samp$Age,            times = dim(ciff.samp)[1]), 
                            life.sg.status_2 = rep(ciff.samp$life.sg.status, times = dim(ciff.samp)[1]), 
                            capture_loc_2    = rep(ciff.samp$capture_loc,    times = dim(ciff.samp)[1]), 
                            year_samp_2      = rep(ciff.samp$year_samp,      times = dim(ciff.samp)[1]), 
                            Month_2          = rep(ciff.samp$Month,          times = dim(ciff.samp)[1]), 
                            year_2           = rep(ciff.samp$year,           times = dim(ciff.samp)[1]), 
                            Number_2         = rep(ciff.samp$Number,         times = dim(ciff.samp)[1]),
                            hsps = 0, 
                            pair_status = paste0(rep(ciff.samp$life.sg.status, each = dim(ciff.samp)[1]),
                                                 rep(ciff.samp$life.sg.status, times = dim(ciff.samp)[1])))
# -------------------
# OK now make the kin
# -------------------
head(ciff.samp.mlt)
ciff.samp.mlt$pair_numbers <- paste(ciff.samp.mlt$Number_1, ciff.samp.mlt$Number_2, sep = "_")
kin_pairs_numbers <- paste(kin.hsps[, 1], kin.hsps[, 2], sep = "_")
check.hsps <- ciff.samp.mlt[which(ciff.samp.mlt$pair_numbers %in% kin_pairs_numbers), ]
ciff.samp.mlt[which(ciff.samp.mlt$pair_numbers %in% kin_pairs_numbers), "hsps"] <- 1
#ciff.samp.mlt[which(ciff.samp.mlt$pair_numbers_1 %in% kin_pairs_numbers_2), "hsps"] <- 1
sum(ciff.samp.mlt$hsps)
# ----------------------
# Subsample a few things
# ----------------------
# Just juveniles-subadult comparisons
ciff.samp.mlt.juv.sa <- ciff.samp.mlt[which(ciff.samp.mlt$pair_status == "JJ"  | 
                                              ciff.samp.mlt$pair_status== "JSA" |
                                              ciff.samp.mlt$pair_status == "SAJ" |
                                              ciff.samp.mlt$pair_status == "SASA"), ]
ciff.samp.mlt.juv.sa.ao           <- ciff.samp.mlt.juv.sa[which(ciff.samp.mlt.juv.sa$year_1 <= ciff.samp.mlt.juv.sa$year_2), ]
ciff.samp.mlt.juv.sa.ao$diff_year <- ciff.samp.mlt.juv.sa.ao$year_2 - ciff.samp.mlt.juv.sa.ao$year_1
# -----------------------
check.hsps <- ciff.samp.mlt.juv.sa[which(ciff.samp.mlt.juv.sa$pair_numbers %in% kin_pairs_numbers), ]
dim(check.hsps)
dim(ciff.samp.mlt.juv.sa)
dim(ciff.samp.mlt.juv.sa.ao)
sum(com.kin$n_comp_HSP)
# Same roost
ciff.samp.mlt.juv.sa.ao$same_roost <- 0
ciff.samp.mlt.juv.sa.ao$same_roost[which(ciff.samp.mlt.juv.sa.ao$capture_loc_1 == ciff.samp.mlt.juv.sa.ao$capture_loc_2)] <- 1
# Try and make Marks table
# ------------------------
ciff.hsp.tab <- table(ciff.samp.mlt.juv.sa.ao$diff_year, ciff.samp.mlt.juv.sa.ao$hsps, ciff.samp.mlt.juv.sa.ao$same_roost)
sum(table(ciff.samp.mlt.juv.sa.ao$diff_year, ciff.samp.mlt.juv.sa.ao$same_roost))
table(ciff.samp.mlt.juv.sa.ao$diff_year, ciff.samp.mlt.juv.sa.ao$pair_status, ciff.samp.mlt.juv.sa.ao$same_roost)
ciff.samp.mlt.juv.sa.ao.hsps <- ciff.samp.mlt.juv.sa.ao[which(ciff.samp.mlt.juv.sa.ao$hsps == 1), ]
table(ciff.samp.mlt.juv.sa.ao.hsps$diff_year, ciff.samp.mlt.juv.sa.ao.hsps$pair_status, ciff.samp.mlt.juv.sa.ao.hsps$same_roost)


sum(ciff.hsp.tab[,,1]) + sum(ciff.hsp.tab[,,2])
ciff.hsp.tab[,,1]
ciff.hsp.tab[,,2]
