# Reusable themes and map elements for making plots for the BIOSHIFTSxFISHERIES manuscript


#' BIOSHIFTS plot theme 1
#' @description 
#' @returns A basic styled theme for all histograms and descriptive plots for 
#' the manuscript
bioshift_plot_theme_1 <- theme_classic() + 
  theme(axis.title=element_text(size=12),
        axis.text = element_text(size=10),
        legend.text=element_text(size=10),
        plot.caption=element_text(size=8),
        plot.title=element_text(size=12, hjust=0.5,),
        plot.subtitle = element_text(size=12, hjust=0.5,),
        strip.text = element_text(size = 12),
        panel.grid.major = element_line(colour="grey", size=0.25),
        strip.background = element_blank())
