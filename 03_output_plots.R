# Run 01_get_protein_info and 02_get_proteoform info first!

# Make ggplots --------------------------------------------------------------------------------

# pb <- progress_bar$new(
#   format = " Creating ggplots [:bar] :percent eta: :eta :spin",
#   total = (length(results_protein) - 2),
#   clear = FALSE, width= 60)

for (i in head(seq_along(results_protein), - 1)) {
  
  #
  #
  
  protein_plot <- ggplot() +
    geom_histogram(data = results_protein[[i]], aes(ave_mass),
                   binwidth = 1000, col = "black", fill = "blue") +
    labs(x = "Average Mass (Da)",
         y = "Protein Count",
         title = basename(names(proteinlist)[i])) +
    xlim(0, 100000) +
    theme_minimal()
  
  ggname <- glue("output/png/{basename(names(proteinlist)[i])}_proteins.png")
  
  ggsave(filename = ggname, plot = protein_plot, 
         device = "png", height = 5.4, width = 9.6, dpi = 600)
  
  
  if (tdreport_file == TRUE) {
    
    
    #
    #
    
    proteoform_plot <- ggplot() +
      geom_histogram(data = results_proteoform[[i]], aes(AverageMass),
                     binwidth = 1000, col = "black", fill = "blue") +
      labs(x = "Average Mass (Da)",
           y = "Proteoform Count",
           title = basename(names(proteoformlist)[i])) +
      xlim(0, 100000) +
      theme_minimal()
    
    ggname <- glue("output/png/{basename(names(proteoformlist)[i])}_proteoforms.png")
    
    ggsave(filename = ggname, plot = proteoform_plot, 
           device = "png", height = 5.4, width = 9.6, dpi = 600)
    
    
    #
    #
    
    ggblap <- ggarrange(protein_plot, proteoform_plot, ncol = 2)
    
    ggname <- glue("output/png/{basename(names(proteinlist)[i])}_all.png")
    
    ggsave(filename = ggname, plot = ggblap,
           device = "png", height = 5.4, width = 9.6, dpi = 600)
    
  }
  
  pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_proteins.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(protein_plot)
  dev.off()
  
  pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_proteoforms.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(proteoform_plot)
  dev.off()
  
  # pb$tick()
  
  }

# Save Workspace ------------------------------------------------------------------------------

# Just in case I want to see an image from a particular run of results

glue("output/workspace_images/{systime}_workspace_image.RData") %>% 
  save.image()

# Testing Other Output Formats ----------------------------------------------------------------
# 
# library(svglite)
# 
# svglite(file = "svglite_test1.svg", width = 9.6, height = 5.4, bg = "transparent", pointsize = 12)
# proteoform_plot
# dev.off()
# 
# library(grDevices)
# 
# pdf(file = "pdf_test1.pdf", width = 9.6, height = 5.4, bg = "transparent")
# proteoform_plot
# dev.off()
# 
# ggsave(filename = "png_test1.png", plot = proteoform_plot, 
#        device = "png", height = 5.4, width = 9.6, dpi = 600)

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
