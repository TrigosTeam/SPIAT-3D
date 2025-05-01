plot_cross_L_gradient3D <- function(cross_L_gradient_df) {
  
  target_cell_types <- colnames(cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient", x = "Radius", y = "Cross L-function value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
