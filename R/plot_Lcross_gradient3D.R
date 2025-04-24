plot_Lcross_gradient3D <- function(Lcross_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(Lcross_gradient_df, "radius", c("observed_Lcross", "expected_Lcross", "Lcross_ratio"))
  plot_result <- plot_result[plot_result$variable != "Lcross_ratio", ]
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Lcross gradient", x = "Radius", y = "Lcross value") +
    scale_colour_discrete(name = "", labels = c("Observed Lcross", "Expected CSR Lcross")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
  
  return(fig) 
}
