calculate_neighbourhood_counts_gradient3D <- function(spe, 
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
    neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radii[i],
                                                                feature_colname,
                                                                FALSE,
                                                                FALSE)
    
    if (is.null(neighbourhood_counts_df)) return(NULL)
    
    neighbourhood_counts_df$ref_cell_id <- NULL
    result[i, ] <- apply(neighbourhood_counts_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_neighbourhood_counts_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
