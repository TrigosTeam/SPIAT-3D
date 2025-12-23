plot_neighbourhood_entropy_gradient3D <- function(neighbourhood_entropy_gradient_df, 
                                                  reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(neighbourhood_entropy_gradient_df, id.vars = c("radius"))
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
    geom_point() +
    geom_line() +
    labs(title = "Neighbourhood entropy gradient", x = "Radius", y = "Neighbourhood entropy", color = "Cell type") +
    theme_bw() +
    ylim(0, 1)
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(neighbourhood_entropy_gradient_df)[seq(ncol(neighbourhood_entropy_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
