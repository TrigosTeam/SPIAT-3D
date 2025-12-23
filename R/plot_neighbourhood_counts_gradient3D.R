plot_neighbourhood_counts_gradient3D <- function(neighbourhood_counts_gradient_df, 
                                                 reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(neighbourhood_counts_gradient_df, "radius")
  
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
    geom_line() + 
    labs(title = "Average neighbourhood counts gradient", x = "Radius", y = "Average neighbourhood counts") + 
    scale_color_discrete(name = "Cell type") +
    theme_bw()
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, sep = ""))
  }
  
  return(fig)
}
