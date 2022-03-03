# =================================================================
# CIFFs simulator function
#   - Simulates a CIFF like population for a specific mortality 
#     rate and saves the population and those who died in 
#     each of the generations to a directory cwd/populations/mort{mratepar}
#   - 40 generations are simulated
#   - A plot of the change in population size is also saved to
#     the directory 
# Author: Luke Lloyd-Jones
# Date started: 04/10/2019
# Date updated: 04/10/2019
# =================================================================
ciff_simulate <- function(start.size, generations, mratepar)
{
  # PARS
  #   start.size  - Initial starting size of population. Needs to be
  #                 thought about when the population is decreasing
  #   generations - how many generations do you want to simulate
  #   mratepar    - mortality rate parameter which stipulates the
  #                 how the population will change
  # RETURNS
  #               - writes population at each generation
  #                 to current working directory
  #               - produces a plot of the total population size
  #                 through time
  # ---------------------------
  # Initial founding population
  # ---------------------------
  # Names go through alphabet for each generation 
  inds.id        <- as.character(paste0("A", seq(1, start.size)))
  # These proportions reflect those that were sampled in reality for CIFFs
  # 2000 is the actual population abundance estimate for CIFFs
  j.size  <- floor(start.size * (1163 / 2000))
  sa.size <- floor(start.size * (279  / 2000))
  a.size  <- (start.size - j.size - sa.size)                       
  life.sg.status <- as.character(rep(c("J", "SA", "A"), 
                                 times = c(j.size, 
                                           sa.size, 
                                           a.size)))
  adult.ages     <- floor(a.size * c(0.7, 0.2, 0.1, 0, 0, 0, 0, 0, 0))
  diff           <- a.size  - sum(adult.ages)  
  ages           <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 
                        times = c(j.size + diff, 
                                  sa.size, 
                                  adult.ages))
  mf             <- rbinom(start.size, 1, 0.5) # 1 for male
  FP             <- "GOD"
  MP             <- "GOD"
  year           <- 1979
  founders       <- data.frame(inds.id, life.sg.status, 
                               ages, mf, FP, MP, year, stringsAsFactors = F)
  population     <- founders
  # ------------------------
  # Cycle over generations
  # ------------------------
  # Make the letters long enough for 40 generations
  LETTERS2        <- c(LETTERS, paste0(LETTERS, LETTERS), paste0(LETTERS, LETTERS, LETTERS))
  pop.samp        <- data.frame(matrix(0, ncol = length(names(population)) + 1))
  names(pop.samp) <- c(names(population), "year_samp")
  years.seq    <- array(0, generations)
  pop.sizes    <- array(0, generations)
  # Start in 1978 so we end up at 2018
  years.seq[1] <- 1978
  pop.sizes[1] <- dim(population)[1]
  # Setup writing directory
  dir.create(file.path(paste0("populations/mort_", mratepar)))
  # Cycle over generations and produce simulated data
  for (gen in seq(2, generations))
  {
    year.i <- 1978 + gen
    years.seq[gen] <- year.i
    print("=========================")
    print(paste0("GENERATION ", gen))
    print("=========================")
    
    # ------------------------------
    # Mate - first year we just add the maximum 209 pairs. Less likely to mate
    # ------------------------------
    adults    <- population[which(population$life.sg.status == "A"), ]
    a.males   <- adults[which(adults$mf == 1), ]
    a.females <- adults[which(adults$mf == 0), ]
    breed.mal     <- a.males[, 1]
    breed.fem     <- a.females[, 1]
    if (length(breed.mal) >  length(breed.fem))
    {
      breed.mal <- sample(breed.mal, length(breed.fem))
    } else if (length(breed.fem) >  length(breed.mal))
    {
      breed.fem <- sample(breed.fem, length(breed.mal))
    }
    # ---
    # Age
    # ---
    population$ages <- population$ages + 1
    # -------------
    # Change status
    # -------------
    population$life.sg.status[which(population$life.sg.status == "SA")] <- "A"
    population$life.sg.status[which(population$life.sg.status == "J")]  <- "SA"
    
    
    # -----------------------------------------------------
    # Recruitment - we have the maximum join the population
    # -----------------------------------------------------
    
    # --------------------
    # Next generation born
    # --------------------
    #kappa          <- 1.1
    #new.nos        <- dim(population)[1] - dim(population)[1] * exp(-kappa)
    inds.id        <- paste0(LETTERS2[gen], "_", seq(1, length(breed.fem)))
    life.sg.status <- array("J", length(inds.id))
    mf             <- rbinom(length(inds.id), 1, 0.5) # 1 for male
    ages           <- array(1, length(inds.id))
    FP             <- breed.fem
    MP             <- breed.mal
    year           <- year.i
    gen.i          <- data.frame(inds.id, life.sg.status, ages, mf, FP, MP, year, stringsAsFactors = F)
    print(paste0("Total new this year ", dim(gen.i)[1]))
    # Bind the new generation with the population
    population <- rbind(gen.i, population)
    hist(population$ages, main = year)
    print(paste0("Population size ", dim(population)[1]))
    pop.sizes[gen] <- dim(population)[1]
    print(paste0("Age distribution "))
    print(table(population$ages))
    print(paste0("Stage distribution "))
    print(table(population$life.sg.status))
    
    # ----------------------------------------------------
    # Die - remove individuals from the current population 
    # ----------------------------------------------------
    # ----------------------------------
    # Given population size kill as per
    # mortality parameter
    # ----------------------------------
    pop.size   <- dim(population)[1]
    alpha      <- mratepar
    dead.nos   <- pop.size - round(pop.size * exp(-alpha), 0)
    dead.samp  <- sample(seq(1, dim(population)[1]), dead.nos)
    dead.samp  <- dead.samp[order(dead.samp)]
    dead       <- population[c(dead.samp), ]
    print(paste0("Total dead this year ", dead.nos , " - RIP"))
    # -----------------------------
    # Write out the dead population
    # -----------------------------
    write.table(dead, paste0("populations/mort_", mratepar, "/ciff_population_gen", gen, "_dead.txt"), 
                col.names = T,
                row.names = F, 
                quote = F, 
                sep = "\t")
    # -----------------
    # Subset population
    # -----------------
    population <- population[-which(population[, 1] %in% as.character(dead[, 1])), ]
    
    write.table(population, paste0("populations/mort_", mratepar, "/ciff_population_gen", gen, ".txt"), 
                col.names = T,
                row.names = F, 
                quote = F, 
                sep = "\t")
  }
# Plot and write out
jpeg(filename = paste0("populations/mort_", mratepar, "/ciff_population_change.jpg"), 
     pointsize = 26, 
     quality = 1000, 
     bg = "white", res = NA, width = 1200, height = 1080)
plot(years.seq, pop.sizes, ylim = c(min(pop.sizes) - 50, max(pop.sizes) + 50))
dev.off() 
}