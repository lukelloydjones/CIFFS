plot_abund_cv <- function(years, N, abund_cv, stderr, lethal, min = 1, max = 1, scale = 1)
{
  require(ggwebthemes)
  # PARS
  #  years     - years to plot over
  #  N         - abundance calculated from pars
  #  abund_cv  - cv of the abundance at each year
  #  stderr    - standard error of estimate on log scale
  #  lethal    - 1 => lethal sampling, 0 => mark-recapture
  #  scale     - Whether abundance is on log scale or natural scale
  # RETURNS
  #  print the plot and saves a version to the results/ subdirectory
  
  #barplot(abund_cv, xlab = "Years", ylab = "CV") # That's cool 
  if (scale == 0)
  {
    abund_res <- data.frame(years, exp(N), 
                            exp(N - 1.96 * stderr), 
                            exp(N + 1.96 * stderr), 
                            exp(N - 1.281552 * stderr), 
                            exp(N + 1.281552 * stderr), 
                            abund_cv)
  } else { 
    abund_res <- data.frame(years, N, 
                            N - 1.96 * stderr, 
                            N + 1.96 * stderr, 
                            N - 1.281552 * stderr, 
                            N + 1.281552 * stderr, 
                            abund_cv)
  }
  colnames(abund_res) <- c("Years", "Abundance", 
                           "Lower", "Upper",
                           "Lower_80", "Upper_80",
                           "CV")
  
  # Set plot axes
  if (min == 1)
  {
    y.ax.min <- min(min(abund_res$Lower), 0) - 50
  } else {
    y.ax.min <- min
  }
  if (max == 1)
  {
    y.ax.max <- min(max(abund_res$Upper), 1e5) + 150
  } else {
    y.ax.max <- max
  }
  
  if (lethal == 1)
  {
    title.plot <- "Lethal sampling"
  } else {
    title.plot <- "Non-lethal sampling"
  }

  # Draw
  a.vals <- ggplot(data = abund_res, aes(x = Years, y = Abundance)) + theme_web_bw() +
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, fill = "grey") +
            geom_ribbon(aes(ymin = Lower_80, ymax = Upper_80), alpha = 0.3, fill = "grey") +
            geom_text(label   = round(abund_res$CV, 3), fontface = 2, size = 5.5,  angle = 0,
                      nudge_y = abund_res$Upper - abund_res$Abundance + 150) +
            geom_line(lwd = 2, colour = "black") +
            xlab("Year")  +
            ylab("Female adult abundance (count)") +
            theme(text  = element_text(size = 25, face = "bold"), 
                  legend.position = "none",
                  axis.text.x  = element_text(size = 15, face = "bold"),
                  axis.title.x = element_text(size = 20, face = "bold"),
                  axis.text.y  = element_text(size = 15, face = "bold"),
                  axis.title.y  = element_text(size = 20, face = "bold")) +
            scale_x_continuous(expand = c(0.025, 0.025), breaks = scales::pretty_breaks(n = dim(abund_res)[1])) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, y.ax.max))
  #ggtitle(label = title.plot, 
  #        subtitle = paste("Over", years.samp, "years")) + 
  print(a.vals)
  
  #brewer.pal(4, "Pastel2")[4]
  
  # Write out
  # jpeg(filename = paste0("results/abundance_ests_plot_scen_", scen, ".jpeg"), pointsize = 26, 
  #      quality = 1000, 
  #      bg = "white", res = NA, width = 1200, height = 1080)
  # print(a.vals)
  # dev.off() 
  return(abund_res)
}
