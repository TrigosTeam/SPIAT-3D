calculate_Gcross_gradient3D <- function(spe, 
                                        reference_cell_type, 
                                        target_cell_type, 
                                        radii, 
                                        feature_colname = "Cell.Type",
                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_Gcross", 
                        "expected_Gcross")
  
  for (i in seq(length(radii))) {
    Gcross_df <- calculate_Gcross3D(spe,
                                    reference_cell_type,
                                    target_cell_type,
                                    radii[i],
                                    feature_colname)
    
    result[i, ] <- Gcross_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_Gcross_gradient3D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
