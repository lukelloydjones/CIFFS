# -----------------------------------------------------
# Function to sample the population of individuals
# from the simulated population
# -----------------------------------------------------
pop_sample <- function(sam.no.py, years.samp, scen, gens, pop.samp.dead, repi)
{
  # PARS
  #  sam.no.py  - number of individuals to sample from each year
  #  years.samp - number of years to sample
  #  scen       - what simulation scenario to sample. Defined by mortality parameter
  #  generations - total number of generations
  #  pop.samp.dead - sample them dead or alive. 1 => dead, 0 => mark recapture
  #  repi - what simulation rep are we doing
  # RETURNS
  #   Writes the sampled populations to sampled_populations/mort_{mrtpar}
  print(paste0("Sampling simulation replicate ", repi))
  # ---------------------------
  # Initial starting population
  # ---------------------------
  # Names go through alphabet for each generation and given number for each individual
  clnames  <- c("inds.id", "life.sg.status", "ages", "mf", "FP", "MP", "year", "year_samp")
  pop.samp <- data.frame(matrix(0, ncol = length(clnames)))
  colnames(pop.samp) <- clnames
  # -----------------------------------------
  # Cycle over years and resample popualation
  # -----------------------------------------
  for (gen in seq(2, gens))
  {
    if (pop.samp.dead) 
    {
      population <- read.table(paste0("populations/mort_", scen, "/ciff_population_gen", gen, "_dead.txt"), 
                               sep = "\t", header = T)
    } else {
      population <- read.table(paste0("populations/mort_", scen, "/ciff_population_gen", gen, ".txt"), 
                               sep = "\t", header = T)
    }
    year          <- 1978 + gen
    if (gen > gens - years.samp)
    {
      print(year)
      samp.no      <- sample(seq(1, dim(population)[1]), sam.no.py)
      samp.no.o    <- samp.no[order(samp.no)]
      pop.samp.gen <- population[samp.no.o, ]
      pop.samp.gen$year_samp <- year
      pop.samp <- rbind(pop.samp, pop.samp.gen)
    }
  }
  # Return the population
  write.table(pop.samp[-1, ], paste0("sampled_populations/mort_", scen, "/pop_samp_", repi, "_yrs_samp_",
                                     years.samp, "_samp_no_py_", sam.no.py, ".txt"),
              col.names = T, row.names = F, quote = F, sep = " ")
}

