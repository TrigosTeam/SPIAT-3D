calculate_Lcross_gradient3D <- function(spe, 
                                        reference_cell_type, 
                                        target_cell_type, 
                                        radii, 
                                        feature_colname = "Cell.Type",
                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 3))
  colnames(result) <- c("observed_Lcross", 
                        "expected_Lcross",
                        "Lcross_ratio")
  
  for (i in seq(length(radii))) {
    Lcross_df <- calculate_Lcross3D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- Lcross_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_Lcross_gradient3D(result, reference_cell_type, target_cell_type)
    fig2 <- plot_Lcross_gradient_ratio3D(result, reference_cell_type, target_cell_type)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}