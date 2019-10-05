# Run 01_get_protein_info and 02_get_proteoform info first!

# Make ggplots --------------------------------------------------------------------------------

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
  
  if (tdreport_file == TRUE) {
    
    pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_proteoforms.pdf"),
        width = 8, height = 5, bg = "transparent")
    print(proteoform_plot)
    dev.off()
    
  }
  
}

# Save Workspace ------------------------------------------------------------------------------

# Just in case I want to see an image from a particular run of results

glue("output/workspace_images/{systime}_workspace_image.RData") %>% 
  save.image()
