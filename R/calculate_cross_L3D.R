calculate_cross_L3D <- function(spe, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K3D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  
  return(result)
}
