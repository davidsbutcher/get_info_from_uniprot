# Run 01_get_protein_info and 02_get_proteoform info first!


# Packages ------------------------------------------------------------------------------------

library(ggpubr)

# Make ggplots --------------------------------------------------------------------------------

for (i in head(seq_along(results_protein), -1)) {
  
  setwd(here("output/png"))
  
  #
  #
  
  protein_plot <- ggplot() +
    geom_histogram(data = results_protein[[i]], aes(ave_mass),
                   binwidth = 1000, col = "black", fill = "blue") +
    labs(x = "Average Mass (Da)",
         y = "Protein Count",
         title = names(proteinlist)[i]) +
    xlim(0, 60000) +
    theme_minimal()
  
  ggname <- glue("{names(proteinlist)[i]}_proteins.png")
  
  ggsave(filename = ggname, plot = protein_plot, device = "png",
         height = 5.4, width = 9.6, dpi = 600)
  
  #
  #
  
  proteoform_plot <- ggplot() +
    geom_histogram(data = results_proteoform[[i]], aes(AverageMass),
                   binwidth = 1000, col = "black", fill = "blue") +
    labs(x = "Average Mass (Da)",
         y = "Proteoform Count",
         title = names(proteoformlist)[i]) +
    xlim(0, 60000) +
    theme_minimal()
  
  ggname <- glue("{names(proteoformlist)[i]}_proteoforms.png")
  
  ggsave(filename = ggname, plot = proteoform_plot, device = "png",
         height = 5.4, width = 9.6, dpi = 600)

  #
  #
  
  ggblap <- ggarrange(protein_plot, proteoform_plot, ncol = 2)
  
  ggname <- glue("{names(proteinlist)[i]}_all.png")
  
  ggsave(filename = ggname, plot = ggblap, device = "png",
         height = 5.4, width = 9.6, dpi = 600)

  setwd(here())
  }

# The number of geom_freqpoly lines cannot exceed the number of 
# input files (until I figure out how to fix it)

ggplot() +
  geom_freqpoly(data = results_proteoform[[1]], aes(AverageMass),
                 binwidth = 1000, col = "black", size = 1.5) +
  geom_freqpoly(data = results_proteoform[[2]], aes(AverageMass),
                binwidth = 1000, col = "blue", size = 1.5) +
  # geom_freqpoly(data = results_proteoform[[3]], aes(AverageMass),
  #               binwidth = 1000, col = "red", size = 1.5) +
  # geom_freqpoly(data = results_proteoform[[4]], aes(AverageMass),
  #               binwidth = 1000, col = "grey", size = 1.5) +
  # geom_freqpoly(data = results_proteoform[[5]], aes(AverageMass),
  #               binwidth = 1000, col = "purple", size = 1.5) +
  labs(x = "Average Mass (Da)",
       y = "Proteoform Count",
       title = "All TD Reports") +
  xlim(0, 60000) +
  theme_minimal()

# Save Workspace ------------------------------------------------------------------------------

# Just in case I want to see an image from a particular run of results

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")

setwd(here("output/workspace_images"))

glue("{systime}_workspace_image.RData") %>% 
  save.image()

setwd(here())


# Deprecated code -----------------------------------------------------------------------------

# temp_plot <- ggplot() +
#   geom_histogram(data = results_protein[[i]], aes(ave_mass),
#                  binwidth = 1000, col = "black", fill = "blue") +
#   # geom_histogram(data = results[[2]], aes(ave_mass),
#   #                binwidth = 1000, col = "black", fill = "red", alpha = 0.4) +
#   # geom_histogram(data = results[[3]], aes(ave_mass),
#   #                binwidth = 1000, col = "black", fill = "yellow", alpha = 0.4) +
#   # geom_freqpoly(data = results_proteoform[[1]], aes(AverageMass), color = "red", binwidth = 1000) +
#   # geom_dotplot(data = results_proteoform[[1]], aes(AverageMass), 
#   #              color = "black", fill = "Red", binwidth = 1000) +
#   # scale_y_continuous() +
#   theme_minimal()
