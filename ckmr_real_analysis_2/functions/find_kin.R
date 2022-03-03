# ----------------------------------------------------
# Function to compute kin in the population sample and 
# write out a melted form of the data
# ----------------------------------------------------
find_kin <- function(sam.no.py, years.samp, scen, repi)
{
  # PARS
  #  sam.no.py  - number of individuals to sample from each year
  #  years.samp - number of years to sample
  #  scen       - what simulation scenario to sample. Defined by mortality parameter
  #  repi - what simulation rep are we doing
  # RETURNS
  #  saves the found pops and hsps for each rep to found_kin/mort_{mrtpar}
  # ---------------------------------------
  # Make a matrix of size n.x.n and compare
  # each individual
  # For each comparison record the relationship if it exists
  # ---------------------------------------
  print(paste0("Finding kin for simulation replicate ", repi))
  pop.samp <- read.table(paste0("sampled_populations/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
                                 years.samp, "_samp_no_py_", sam.no.py, ".txt"), 
                         header = T, stringsAsFactors = F)
  n <- dim(pop.samp)[1]
  pop.samp.mlt <- data.frame(rep(pop.samp$inds.id , each = n), # Parents start
                             rep(pop.samp$life.sg.status , each = n), #
                             rep(pop.samp$ages , each = n), #
                             rep(pop.samp$mf , each = n), #
                             rep(pop.samp$FP , each = n), #
                             rep(pop.samp$MP , each = n), #
                             rep(pop.samp$year , each = n), 
                             rep(pop.samp$year_samp , each = n),# Parents end
                             rep(pop.samp$inds.id, times = n), # Juveniles start
                             rep(pop.samp$life.sg.status, times = n), #
                             rep(pop.samp$ages, times = n), #
                             rep(pop.samp$mf, times = n), #
                             rep(pop.samp$FP, times = n), #
                             rep(pop.samp$MP, times = n), #
                             rep(pop.samp$year, times = n), 
                             rep(pop.samp$year_samp , times = n),# Juveniles end
                             0, # MOP
                             0, # FOP
                             0, # MHSP
                             0, # FHSP
                             stringsAsFactors = F) 
  indices <- data.frame(rep(seq(1, n), each = n),  rep(seq(1, n) , times = n))
  colnames(pop.samp.mlt) <- c(paste0(colnames(pop.samp), "_P"), paste0(colnames(pop.samp), "_O"), "MOP", "FOP", "MHSP", "FHSP")
  print("Computing POPS and HSPS")
  # ----
  # POPS
  # ----
  
  check1 <- which(pop.samp.mlt$inds.id_P == pop.samp.mlt$FP_O) # Is indivdual i j's mother?
  check2 <- which(pop.samp.mlt$inds.id_P == pop.samp.mlt$MP_O) # Is indivdual i j's father?
  
  pop.samp.mlt[check1, "MOP"] = 1
  pop.samp.mlt[check2, "FOP"] = 1
  
  # Add indices to melted data frame
  colnames(indices) <- c("index_parent", "index_child")
  pop.samp.mlt      <- data.frame(pop.samp.mlt, indices)
  
  # Subset to potential mops and fops
  mops.mlt          <- pop.samp.mlt[check1, ]
  fops.mlt          <- pop.samp.mlt[check2, ]
  
  # -------------------------------------------------------------------------------
  # MOPS - This next process should be fine for lethal sampling as it NEVER HAPPENS
  # -------------------------------------------------------------------------------
  # Remove those that are adult-adult mops
  mops.pair.acs <- mops.mlt[-which(paste0(mops.mlt$life.sg.status_P, mops.mlt$life.sg.status_O) == "AA"), ]
  # For mark-recapture scenarion adult second is no good. Remove.
  if (length(which(mops.pair.acs$life.sg.status_O == "A")) > 0)
  {
    mops.pair.acs <- mops.pair.acs[-which(mops.pair.acs$life.sg.status_O == "A"), ]
  }
  # Make a data frame just of the pair-ids
  mops.mlt.df      <- data.frame(mops.pair.acs[, "inds.id_P"], mops.pair.acs[, "inds.id_O"], stringsAsFactors = F)
  # Get the unique set of mops from the potential set of duplicates 
  mops.mlt.df.dup  <- unique(mops.mlt.df[which(duplicated(mops.mlt.df)), ])
  mops.mlt.df      <- mops.mlt.df[order(mops.mlt.df[, 1]), ]
  colnames(mops.mlt.df.dup) <- c("inds.id_P", "inds.id_O")
  mops.all.dup     <- right_join(mops.pair.acs, mops.mlt.df.dup, by = c("inds.id_P", "inds.id_O"))
  # Preserve those that were mops but not part of the duplicate madness
  mops.pair.ajs    <- anti_join(mops.pair.acs, mops.all.dup, by = c("inds.id_P", "inds.id_O"))
  # # Take 1st POP
  mops.all.dup.u1      <- mops.all.dup[-which(duplicated(data.frame(mops.all.dup[, "inds.id_P"], mops.all.dup[, "inds.id_O"], stringsAsFactors = F))), ]
  mops.all.dup.u1$diff <- mops.all.dup.u1$year_samp_P - mops.all.dup.u1$year_O
  # Take 2nd POP
  mops.all.dup.u2      <- mops.all.dup[which(duplicated(data.frame(mops.all.dup[, "inds.id_P"], mops.all.dup[, "inds.id_O"], stringsAsFactors = F))), ]
  mops.all.dup.u2$diff <- mops.all.dup.u2$year_samp_P - mops.all.dup.u2$year_O
  mops.all.dup.u2      <- mops.all.dup.u2 %>% distinct(inds.id_P, inds.id_O, .keep_all = T)
  # Do a random sample
  # samp.no          <- sample(seq(1, dim(mops.all.dup)[1]), dim(mops.mlt.df.dup)[1])
  # samp.no.o        <- samp.no[order(samp.no)]
  # mops.all.dup.u   <- mops.all.dup[samp.no.o, ]
  # Take the indices for later use
  #mops.inds        <- rbind(mops.pair.ajs[, c("index_parent", "index_child")], mops.all.dup.u2[, c("index_parent", "index_child")]) 
  # Take all of them (No adjustment)
  mops.inds        <- rbind(mops.pair.ajs[, c("index_parent", "index_child")], mops.all.dup[, c("index_parent", "index_child")])  
  # Take NONE of them
  #mops.inds        <- rbind(mops.pair.ajs[, c("index_parent", "index_child")])  
  
  # -------------------------------------------------------------------------------
  # FOPS - This next process should be fine for lethal sampling as it NEVER HAPPENS
  # -------------------------------------------------------------------------------
  # Remove those that are adult-adult mops
  fops.pair.acs <- fops.mlt[-which(paste0(fops.mlt$life.sg.status_P, fops.mlt$life.sg.status_O) == "AA"), ]
  # For mark-recapture scenarion adult second is no good. Remove.
  if (length(which(fops.pair.acs$life.sg.status_O == "A")) > 0)
  {
    fops.pair.acs <- fops.pair.acs[-which(fops.pair.acs$life.sg.status_O == "A"), ]
  }
  # Make a data frame just of the pair-ids
  fops.mlt.df      <- data.frame(fops.pair.acs[, "inds.id_P"], fops.pair.acs[, "inds.id_O"], stringsAsFactors = F)
  # Get the unique set of fops from the potential set of duplicates 
  fops.mlt.df.dup  <- unique(fops.mlt.df[which(duplicated(fops.mlt.df)), ])
  fops.mlt.df      <- fops.mlt.df[order(fops.mlt.df[, 1]), ]
  colnames(fops.mlt.df.dup) <- c("inds.id_P", "inds.id_O")
  fops.all.dup     <- right_join(fops.pair.acs, fops.mlt.df.dup, by = c("inds.id_P", "inds.id_O"))
  # Preserve those that were fops but not part of the duplicate madness
  fops.pair.ajs    <- anti_join(fops.pair.acs, fops.all.dup, by = c("inds.id_P", "inds.id_O"))
  # Take 1st POP
  #fops.all.dup.u   <- fops.all.dup[-which(duplicated(data.frame(fops.all.dup[, "inds.id_P"], fops.all.dup[, "inds.id_O"], stringsAsFactors = F))), ]
  # Take 2nd POP. Why in god's name does this work (09/12/2019)
  fops.all.dup.u2   <- fops.all.dup[which(duplicated(data.frame(fops.all.dup[, "inds.id_P"], fops.all.dup[, "inds.id_O"], stringsAsFactors = F))), ]
  fops.all.dup.u2   <- fops.all.dup.u2 %>% distinct(inds.id_P, inds.id_O, .keep_all = T)
  # Do a random sample
  # samp.no          <- sample(seq(1, dim(fops.all.dup)[1]), dim(fops.mlt.df.dup)[1])
  # samp.no.o        <- samp.no[order(samp.no)]
  #fops.all.dup.u   <- fops.all.dup[samp.no.o, ]
  # Take the indices for later use
  #fops.inds        <- rbind(fops.pair.ajs[, c("index_parent", "index_child")], fops.all.dup.u2[, c("index_parent", "index_child")]) 
  # Take the indices for later use
  #fops.inds        <- rbind(fops.pair.ajs[, c("index_parent", "index_child")], fops.all.dup.u[, c("index_parent", "index_child")])  
  # Take ALL of them (No adjustment)
  fops.inds        <- rbind(fops.pair.ajs[, c("index_parent", "index_child")], fops.all.dup[, c("index_parent", "index_child")])  
  # Take NONE of them
  #fops.inds        <- rbind(fops.pair.ajs[, c("index_parent", "index_child")])  
  
  # ----
  # HSPS
  # ----
  check3 <- ((pop.samp.mlt$FP_P == pop.samp.mlt$FP_O) & (pop.samp.mlt$inds.id_P != pop.samp.mlt$inds.id_O)) # Is indivdual i j's mother?
  check4 <- ((pop.samp.mlt$MP_P == pop.samp.mlt$MP_O) & (pop.samp.mlt$inds.id_P != pop.samp.mlt$inds.id_O)) # Is indivdual i j's father?
  
  pop.samp.mlt[check3, "MHSP"] = 1
  pop.samp.mlt[check4, "FHSP"] = 1
  
  mhsps.inds <- indices[check3, ]
  phsps.inds <- indices[check4, ]
  # -------------------
  # Tally up some stuff
  # -------------------
  # POPs
  sum(pop.samp.mlt$MOP) #
  sum(pop.samp.mlt$FOP) #
  pops.tot <- sum(pop.samp.mlt$FOP) + sum(pop.samp.mlt$MOP) #Get the same number of POPs and stuff. That's good
  # HSPs
  sum(pop.samp.mlt$MHSP) #
  sum(pop.samp.mlt$FHSP) #
  hsps.tot <- sum(pop.samp.mlt$MHSP) + sum(pop.samp.mlt$FHSP) #Get the same number of POPs and stuff. That's good
  
  kin <- list(mops.inds, fops.inds, mhsps.inds, phsps.inds)
  names(kin) <- c("mops", "fops", "mhsps", "phsps")
  save(kin, file = paste0("found_kin/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
                          years.samp, "_samp_no_py_", sam.no.py, ".Rdata"))
}
