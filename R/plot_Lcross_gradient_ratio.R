plot_Lcross_gradient_ratio3D <- function(Lcross_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- data.frame(radius = Lcross_gradient_df$radius,
                            observed_Lcross_gradient_ratio = Lcross_gradient_df$Lcross_ratio,
                            expected_Lcross_gradient_ratio = 1)
  
  plot_result <- reshape2::melt(plot_result, "radius", c("observed_Lcross_gradient_ratio", "expected_Lcross_gradient_ratio"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Lcross ratio gradient", x = "Radius", y = "Lcross ratio") +
    scale_colour_discrete(name = "", labels = c("Observed Lcross ratio", "Expected CSR Lcross ratio")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
  
  return(fig) 
}
