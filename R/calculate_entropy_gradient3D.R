calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    entropy_df <- calculate_entropy3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    if (is.null(entropy_df)) return(NULL)
    
    result[i, ] <- apply(entropy_df[ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    fig <- plot_entropy_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
