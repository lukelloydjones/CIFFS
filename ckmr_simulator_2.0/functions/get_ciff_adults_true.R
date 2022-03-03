get_ciff_adults <- function(sam.no.py, years.samp, gens, scen, years.tot, two_or_none)
{
  # Get the global truth from the actual populations
  adults.ts <- matrix(0, nrow = years.samp + 2, ncol = 3)
  i=0
  for (gen in seq((gens - years.tot + 1), gens))
  {
    i = i + 1
    print(gen)
    population <- read.table(paste0("populations", two_or_none, "/mort_", scen, "/ciff_population_gen", gen, ".txt"), sep = "\t", header = T)
    population.adults <- population[which(population$life.sg.status == "A"), ]
    adults.ts[i, 1] <- table(population$life.sg.status, population$mf)[1, 1]
    adults.ts[i, 2] <- table(population$life.sg.status, population$mf)[1, 2]
    adults.ts[i, 3] <- adults.ts[i, 1] + adults.ts[i, 2]
  }
  return(list(adults.ts[, 3], adults.ts[, 1], adults.ts[, 2]))
}