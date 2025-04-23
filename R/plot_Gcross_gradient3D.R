plot_Gcross_gradient3D <- function(Gcross_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(Gcross_gradient_df, "radius", c("observed_Gcross", "expected_Gcross"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Gcross gradient", x = "Radius", y = "Gcross value") +
    scale_colour_discrete(name = "", labels = c("Observed Gcross", "Expected CSR Gcross")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
  
  return(fig) 
}
