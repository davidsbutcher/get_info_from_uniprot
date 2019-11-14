# Run 01_get_protein_info and 02_get_proteoform info first!

if (dir.exists("output/png") == FALSE) dir.create("output/png")

if (dir.exists("output/pdf") == FALSE) dir.create("output/pdf")

# Make ggplots --------------------------------------------------------------------------------

for (i in head(seq_along(results_protein), - 1)) {
  
  ##### Mass histograms for protein results, 1000 m/z bins
  
  protein_plot <- 
    ggplot() +
    geom_histogram(data = results_protein[[i]], aes(ObservedPrecursorMass),
                   binwidth = 1000, col = "black", fill = "blue") +
    labs(x = "Average Mass (Da)",
         y = "Protein Count",
         title = basename(names(proteinlist)[i])) +
    xlim(0,
         plyr::round_any(
           max(results_protein[[i]]$ObservedPrecursorMass),
           10000, f = ceiling)
    ) +
    theme_minimal()
  
  ##### Save protein mass histograms
  
  ggname <- glue("output/png/{basename(names(proteinlist)[i])}_proteins.png")
  
  ggsave(filename = ggname, plot = protein_plot, 
         device = "png", height = 5.4, width = 9.6, dpi = 600)
  
  pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_proteins.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(protein_plot)
  dev.off()
  
  ##### Heatmaps for protein results, Best Hits
  
  bin_size <- 1000   # SPECIFY SIZE for heatmap mass bins
  bins <-  seq.int(0, 200000, by = bin_size)
  
  protein_heatmap_data <- 
    results_protein[[i]] %>% 
    select(c("ObservedPrecursorMass", "fraction")) %>% 
    mutate(mass_bin = cut(ObservedPrecursorMass, breaks = bins, labels = FALSE)) %>%
    group_by(fraction, mass_bin) %>% 
    summarize(n()) %>%
    mutate(mass_bin =  (mass_bin * bin_size) - bin_size) %>% 
    rename("Protein Count" = `n()`) %>%
    ungroup() %>% 
    mutate(fraction = forcats::as_factor(fraction))
  
  protein_heatmap <- 
    protein_heatmap_data %>% 
    ggplot(aes(mass_bin, fraction,
               fill = `Protein Count`)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", direction = -1) +
    xlim(0,
         plyr::round_any(
           max(protein_heatmap_data$mass_bin),
           10000, f = ceiling)
    ) +
    scale_y_discrete(limits = rev(levels(protein_heatmap_data$fraction))) +
    labs(x = "Mass Bin (Da)",
         y = "Fraction") +
    theme_minimal() +
    theme(text = element_text(size=18))
  
  ##### Save protein heatmaps, Best Hits
  
  heatmapname <- 
    glue("output/png/{basename(names(proteinlist)[i])}_protein_heatmap.png")
  
  ggsave(filename = heatmapname, plot = protein_heatmap, 
         device = "png", height = 5.4, width = 9.6, dpi = 600)
  
  pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_protein_heatmap.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(protein_heatmap)
  dev.off()
  
  ##### Heatmaps for protein results, ALL HITS
  
  bin_size <- 1000   # SPECIFY SIZE for heatmap mass bins
  bins <-  seq.int(0, 200000, by = bin_size)
  
  protein_heatmap_allhits_data <- 
    proteinlistfull[[i]] %>% 
    select(c("ObservedPrecursorMass", "fraction")) %>% 
    mutate(mass_bin = cut(ObservedPrecursorMass, breaks = bins, labels = FALSE)) %>%
    group_by(fraction, mass_bin) %>% 
    summarize(n()) %>%
    mutate(mass_bin =  (mass_bin * bin_size) - bin_size) %>% 
    rename("Protein Count" = `n()`) %>%
    ungroup() %>% 
    mutate(fraction = forcats::as_factor(fraction))
  
  protein_heatmap_allhits <- 
    protein_heatmap_allhits_data %>% 
    ggplot(aes(mass_bin, fraction,
               fill = `Protein Count`)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", direction = -1) +
    xlim(0,
         plyr::round_any(
           max(protein_heatmap_allhits_data$mass_bin),
           10000, f = ceiling)
    ) +
    scale_y_discrete(limits = rev(levels(protein_heatmap_allhits_data$fraction))) +
    labs(x = "Mass Bin (Da)",
         y = "Fraction") +
    theme_minimal() +
    theme(text = element_text(size=18))
  
  ##### Save protein heatmaps, ALL Hits
  
  heatmapname <- 
    glue("output/png/{basename(names(proteinlist)[i])}_protein_heatmap_allhits.png")
  
  ggsave(filename = heatmapname, plot = protein_heatmap_allhits, 
         device = "png", height = 5.4, width = 9.6, dpi = 600)
  
  pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_protein_heatmap_allhits.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(protein_heatmap_allhits)
  dev.off()
  
  if (tdreport_file == TRUE) {
    
    ##### Mass histograms for proteoform results, 1000 m/z bins 
    
    proteoform_plot <- ggplot() +
      geom_histogram(data = results_proteoform[[i]], aes(AverageMass),
                     binwidth = 1000, col = "black", fill = "blue") +
      labs(x = "Average Mass (Da)",
           y = "Proteoform Count",
           title = basename(names(proteoformlist)[i])) +
      xlim(0, 100000) +
      theme_minimal()
    
    ##### Save proteoform mass histograms
    
    ggname <- glue("output/png/{basename(names(proteoformlist)[i])}_proteoforms.png")
    
    ggsave(filename = ggname, plot = proteoform_plot, 
           device = "png", height = 5.4, width = 9.6, dpi = 600)
    
    pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_proteoforms.pdf"),
        width = 8, height = 5, bg = "transparent")
    print(proteoform_plot)
    dev.off()
    
    ##### Save combined protein and proteoform plots
    
    ggblap <- ggarrange(protein_plot, proteoform_plot, ncol = 2)
    
    ggname <- glue("output/png/{basename(names(proteinlist)[i])}_all.png")
    
    ggsave(filename = ggname, plot = ggblap,
           device = "png", height = 5.4, width = 9.6, dpi = 600)
    
    ##### Heatmaps for proteoform results
    
    bin_size <- 1000   # SPECIFY SIZE for heatmap mass bins
    bins <-  seq.int(0, 200000, by = bin_size)
    
    proteoform_heatmap_data <- 
      results_proteoform[[i]] %>% 
      select(c("ObservedPrecursorMass", "fraction")) %>% 
      mutate(mass_bin = cut(ObservedPrecursorMass, breaks = bins, labels = FALSE)) %>%
      group_by(fraction, mass_bin) %>% 
      summarize(n()) %>%
      mutate(mass_bin =  (mass_bin * bin_size) - bin_size) %>% 
      rename("Proteoform Count" = `n()`) %>%
      ungroup() %>% 
      mutate(fraction = forcats::as_factor(fraction))
    
    proteoform_heatmap <- 
      proteoform_heatmap_data %>% 
      ggplot(aes(mass_bin, fraction,
                 fill = `Proteoform Count`)) +
      geom_tile() +
      scale_fill_viridis_c(option = "C", direction = -1) +
      xlim(0,
           plyr::round_any(
             max(proteoform_heatmap_data$mass_bin),
             10000, f = ceiling)
      ) +
      scale_y_discrete(limits = rev(levels(proteoform_heatmap_data$fraction))) +
      labs(x = "Mass Bin (Da)",
           y = "Fraction") +
      theme_minimal() +
      theme(text = element_text(size=18))
    
    ##### Save proteoform heatmaps
    
    heatmapname <- 
      glue("output/png/{basename(names(proteinlist)[i])}_proteoform_heatmap.png")
    
    ggsave(filename = heatmapname, plot = proteoform_heatmap, 
           device = "png", height = 5.4, width = 9.6, dpi = 600)
    
    pdf(file = glue("output/pdf/{basename(names(proteinlist)[i])}_proteoform_heatmap.pdf"),
        width = 8, height = 5, bg = "transparent")
    print(proteoform_heatmap)
    dev.off()
    
  }
  
}

# Save Workspace ------------------------------------------------------------------------------

# Just in case you want to see an image from a particular run of results

if (dir.exists("output/workspace_images") == FALSE) dir.create("output/workspace_images")

glue("output/workspace_images/{systime}_workspace_image.RData") %>% 
  save.image()
