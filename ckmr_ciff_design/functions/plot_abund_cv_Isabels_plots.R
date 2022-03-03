plot_abund_cv <- function(years, N, abund_cv, lethal, min = 1, max = 1)
{
  # PARS
  #  years     - years to plot over
  #  N         - abundance calculated from pars
  #  abund_cv  - cv of the abundance at each year
  #  lethal    - 1 => lethal sampling, 0 => mark-recapture
  # RETURNS
  #  print the plot and saves a version to the results/ subdirectory
  
  #barplot(abund_cv, xlab = "Years", ylab = "CV") # That's cool 
  
  abund_res <- data.frame(years, N, N - 1.96 * abund_cv * N, N + 1.96 * abund_cv * N, abund_cv)
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
  a.vals <- ggplot(data = abund_res, aes(x = Years, y = Abundance)) + 
            theme_classic() +
            geom_ribbon(aes(ymin = Lower, ymax = Upper, alpha = 0.1, fill = "colours")) +
            geom_text(label = round(abund_res$CV, 3), aes(fontface=2), angle = 35,
                      nudge_y = abund_res$Upper - abund_res$Abundance + 250) +
            geom_line(lwd = 2, colour = "midnightblue") +
            ylim(c(y.ax.min, y.ax.max)) +
            xlab("Year") +
            ylab("Adult abundance estimate") +
            ggtitle(label = "Coefficient of Variation Estimates", #
            subtitle = paste("3 plates, 288 samples")) + 
            scale_fill_manual(values="plum", name="colours") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = dim(abund_res)[1])) +
            theme(text  = element_text(size = 15), legend.position = "none",
                  axis.text  = element_text(colour = "black"),
                  axis.title = element_text(size = 15),
                  axis.title.x  = element_text(vjust=-2),
                  axis.title.y  = element_text(vjust=3),
                  plot.title = element_text(size = 20, face = "bold", vjust=3, hjust=0.5),
                  plot.subtitle = element_text(vjust=2, hjust=0.5),
                  plot.margin=unit(c(20,15,15,20),"pt"),
                  panel.grid.major = element_line(color="grey", linetype="dashed"))
  print(a.vals)
  
  #brewer.pal(4, "Pastel2")[4]
  
  # Write out
  # jpeg(filename = paste0("results/abundance_ests_plot_scen_", scen, ".jpeg"), pointsize = 26, 
  #      quality = 1000, 
  #      bg = "white", res = NA, width = 1200, height = 1080)
  # print(a.vals)
  # dev.off() 
}
