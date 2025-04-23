calculate_Gcross3D <- function(spe,
                               reference_cell_type,
                               target_cell_type,
                               radius,
                               feature_colname = "Cell.Type") {
  
  # Get the number of target cells in the radius around each reference cell
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_type,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  reference_target_interactions <- cells_in_neighbourhood_df[[target_cell_type]]
  
  # Gcross: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  result <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  return(result)
}
