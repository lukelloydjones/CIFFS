plot_abundance <- function(abundance.ests, sam.no.py, years.samp, gens, scen, lethal, min = 1, max = 1)
{
  # PARS
  #  abundance.ests - abundance estimate from across the simulations
  #  sam.no.py  - number of individuals to sample from each year
  #  years.samp - number of years to sample
  #  scen       - what simulation scenario to sample. Defined by mortality parameter
  #  repi - what simulation rep are we doing
  # RETURNS
  #  print the plot and saves a version to the results/ subdirectory
  library(ggplot2)
  source("functions/get_ciff_adults_true.R")
  # Various statistics
  quants  <- apply(abundance.ests, 2, function(x) {quantile(x, probs = c(0.025, 0.975), na.rm = T)})
  means   <- colMeans(abundance.ests, na.rm = T)
  medians <- apply(abundance.ests, 2, median, na.rm = T)
  truth   <- get_ciff_adults(sam.no.py, years.samp, gens, scen, dim(abundance.ests)[2])[[1]]
  # Bind up
  adult.ts.df <- as.data.frame(t((rbind(abundance.ests, quants, means, medians, truth))))
  adult.ts.df$years <- colnames(abundance.ests)
  # Set plot axes
  if (min == 1)
  {
    y.ax.min <- max(min(abundance.ests), 0) - 50
  } else {
    y.ax.min <- min
  }
  if (max == 1)
  {
    y.ax.max <- min(max(abundance.ests), 1e5) + 50
  } else {
    y.ax.max <- max
  }
  
  if (lethal == 1)
  {
    title.plot <- "Lethal sampling"
  } else {
    title.plot <- "Non-lethal sampling"
  }
  # Melt for ggplot
  adult.ts.df.m <- melt(adult.ts.df, id = "years")
  adult.ts.df.m[, 2] <- as.character(adult.ts.df.m[, 2])
  adult.ts.df.m$Extra <- adult.ts.df.m[, 2]
  adult.ts.df.m[grep("V", adult.ts.df.m[, "Extra"]), "Extra"] <- "SIM"
  adult.ts.df.m[grep("%", adult.ts.df.m[, "Extra"]), "Extra"] <- "QUANT"
  adult.ts.df.m$Extra <- factor(adult.ts.df.m$Extra)
  # Reorder levels for colour specificity
  levels(adult.ts.df.m$Extra) <- c("Mean", "Median", "Quantiles", "Replicates", "Truth")#, "truth2")
  adult.ts.df.m$Extra <- factor(adult.ts.df.m$Extra, levels = c("Replicates", "Mean", "Median", "Quantiles", "Truth"))#, "truth2"))
  adult.ts.df.m$variable <- factor(adult.ts.df.m$variable, levels = c(paste0("V", seq(1,30)), "2.5%", 
                                                                   "97.5%", 
                                                                   "means", "medians", "truth"))
  # Draw
  a.vals <-  ggplot(data = adult.ts.df.m, aes(x = years, y = value, group = variable)) + 
    theme_classic() +
    geom_line(lwd = 3, aes(color = Extra, alpha = Extra)) +
    ylim(c(y.ax.min, y.ax.max)) +
    xlab("Year") +
    ylab("Adult abundance estimate") +
    ggtitle(label = "CMKR Simulator Results",
            subtitle = paste0(title.plot, "; ", sam.no.py, " samples per year, over ", years.samp, " years")) + 
    scale_alpha_manual(values = c(0.15, 1, 1, 1, 1, 1), name = "Extra") +
    scale_color_manual(values = c(brewer.pal(7, "Set3")[c(5, 4, 7, 3, 6)]), name = "Extra") +
    theme(text  = element_text(size = 35), legend.title = element_blank(),
          axis.text  = element_text(colour = "black"),
          axis.title = element_text(size = 35),
          axis.title.x  = element_text(vjust=-2),
          axis.title.y  = element_text(vjust=3),
          plot.title = element_text(size = 40, face = "bold", vjust=3, hjust=0.5),
          plot.subtitle = element_text(vjust=4),
          plot.margin=unit(c(50,40,40,50),"pt"),
    panel.grid.major = element_line(color="grey", linetype="dashed"))
  print(a.vals)
  # Write out
  # jpeg(filename = paste0("results/isabels_plot_scen_", scen, ".jpeg"), pointsize = 26, 
  #      quality = 1000, bg = "white", res = NA, width = 1200, height = 800)
  png(filename = paste0("results/isabels_plot_scen_", scen, ".png"), pointsize = 26, 
       bg = "white", res = NA, width = 1400, height = 900)
  print(a.vals)
  dev.off() 
}
