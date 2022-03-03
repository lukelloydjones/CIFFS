plot_abund_cv <- function(years, N, abund_cv, std.err, lethal, min = 1, max = 1)
{
  # PARS
  #  years     - years to plot over
  #  N         - abundance calculated from pars
  #  abund_cv  - cv of the abundance at each year
  #  std.err   - standard error of abundance estimate
  #  lethal    - 1 => lethal sampling, 0 => mark-recapture
  # RETURNS
  #  print the plot and saves a version to the results/ subdirectory
  
  #barplot(abund_cv, xlab = "Years", ylab = "CV") # That's cool 
  
  abund_res <- data.frame(years, exp(N), exp(N - 1.96 * std.err), exp(N + 1.96 * std.err), abund_cv)
  colnames(abund_res) <- c("Years", "Abundance", "Lower", "Upper", "CV")
  
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
  a.vals <- ggplot(data = abund_res, aes(x = Years, y = Abundance)) + theme_grey() +
            geom_ribbon(aes(ymin = Lower, ymax = Upper, alpha = 0.1, fill = "colours")) +
            geom_text(label = round(abund_res$CV, 3), aes(fontface=2), size = 5, angle = 35,
                      nudge_y = abund_res$Upper - abund_res$Abundance + 150) +
            geom_line(lwd = 2, colour = "midnightblue") +
            ylim(c(y.ax.min, y.ax.max)) + ylab("") + xlab("") + 
            theme(text  = element_text(size = 25, face = "bold"), 
                  legend.position = "none",
                  axis.text.x  = element_text(size = 15, face = "bold"),
                  axis.text.y  = element_text(size = 15, face = "bold")) +
            scale_fill_manual(values="plum", name="colours") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = dim(abund_res)[1])) 
  print(a.vals)
  #            ggtitle(label = title.plot, 
  # subtitle = paste("Over", years.samp, "years")) + 
  # brewer.pal(4, "Pastel2")[4]
  # xlab("Years modelled") +
  # ylab("Adult abundance (count)") +
  
  # Write out
  # jpeg(filename = paste0("results/abundance_ests_plot_scen_", scen, ".jpeg"), pointsize = 26, 
  #      quality = 1000, 
  #      bg = "white", res = NA, width = 1200, height = 1080)
  # print(a.vals)
  # dev.off() 
}
