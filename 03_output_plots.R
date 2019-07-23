# Run 01_get_protein_info and 02_get_proteoform info first!

results_plot <- ggplot() +
  # geom_histogram(data = results[[1]], aes(ave_mass),
  #                binwidth = 1000, col = "black", fill = "blue", alpha = 0.4) +
  # geom_histogram(data = results[[2]], aes(ave_mass),
  #                binwidth = 1000, col = "black", fill = "red", alpha = 0.4) +
  # geom_histogram(data = results[[3]], aes(ave_mass),
  #                binwidth = 1000, col = "black", fill = "yellow", alpha = 0.4) +
  geom_freqpoly(data = results_proteoform[[1]], aes(ave_mass), color = "red", binwidth = 1000) +
  theme_minimal()