# Cases to test: ------
# 1. Correct spe input
md1 <- spe_metadata_background_template("random")
md1 <- spe_metadata_cluster_template("regular", "sphere", md1)
md1 <- spe_metadata_cluster_template("regular", "ellipsoid", md1)
md1 <- spe_metadata_cluster_template("regular", "cylinder", md1)
spe1 <- simulate_spe_metadata3D(md1)
spe1$Cell.Type[spe1$Cell.Type == "Endothelial"] <- "C"
spe1$Cell.Type[spe1$Cell.Type == "Tumour"] <- "A"
spe1$Cell.Type[spe1$Cell.Type == "Immune"] <- "B"
spe1$Cell.Type[spe1$Cell.Type == "Others"] <- "O"

# 2. Not spe object input
spe2 <- list()

# 3. Spe object input, but with invalid feature colname
md3 <- spe_metadata_background_template("random")
spe3 <- simulate_spe_metadata3D(md3)
names(colData(spe3))[1] <- "dummy"

# 4. Zero cells
spe4 <- SpatialExperiment()

# 5. One cell
md5 <- spe_metadata_background_template("random")
md5$background$n_cells <- 1
md5$background$cell_types <- c("A", "B")
spe5 <- simulate_spe_metadata3D(md5)

# 6. Two cells
md6 <- spe_metadata_background_template("random")
md6$background$n_cells <- 2
md6$background$cell_types <- c("A", "B")
spe6 <- simulate_spe_metadata3D(md6)

# 7. Many cells, but 0 cell type of interest/reference/target
md7 <- spe_metadata_background_template("random")
md7$background$cell_types <- c("C")
md7$background$cell_proportions <- c(1)
spe7 <- simulate_spe_metadata3D(md7)

# 8. Many cells, but 1 cell type of interest/reference/target
md8 <- spe_metadata_background_template("random")
md8$background$cell_types <- c("C")
md8$background$cell_proportions <- c(1)
spe8 <- simulate_spe_metadata3D(md8)
spe8[["Cell.Type"]][1] <- "A"

# 9. Many cells, but 2 cell type of interest/reference/target
md9 <- spe_metadata_background_template("random")
md9$background$cell_types <- c("C")
md9$background$cell_proportions <- c(1)
spe9 <- simulate_spe_metadata3D(md9)
spe9[["Cell.Type"]][1] <- "A"
spe9[["Cell.Type"]][1] <- "B"

# 10. Reference cell type and target cell type are separated.
md10 <- spe_metadata_background_template("random")
md10$background$cell_types <- c("C")
md10$background$cell_proportions <- c(1)
spe10 <- simulate_spe_metadata3D(md10)
spe10$Cell.Type[spatialCoords(spe10[ , 3]) < 20] <- "A"
spe10$Cell.Type[spatialCoords(spe10[ , 3]) > 280] <- "B"

# Put all spes into a list
spes <- list(spe1,
             spe2,
             spe3,
             spe4,
             spe5,
             spe6,
             spe7,
             spe8,
             spe9,
             spe10)

# Calculate_cell_proportions3D -------

error_df <- data.frame(matrix(nrow = 400, ncol = 5))
colnames(error_df) <- c("spe", "cell_type_of_interest", "feature_colname", "plot_image", "error_message")
index <- 1

cell_types_of_interest_options <- list(NULL, "A", c("A", "B"), c("A", "B", "C"), c("A", "B", "C", "D"), "Invalid")
feature_colname_options <- c("Cell.Type", "Invalid")
plot_image_options <- list(TRUE, FALSE, "Invalid")

for (i in seq(length(spes))) {
  
  curr_spe <- spes[[i]]
  
  for (cell_types_of_interest in cell_types_of_interest_options) {
    error_output <- tryCatch(
      expr = {
        calculate_cell_proportions3D(curr_spe,
                                     cell_types_of_interest = cell_types_of_interest,
                                     feature_colname = "Cell.Type",
                                     plot_image = TRUE)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), ifelse(is.null(cell_types_of_interest), "NULL", paste(cell_types_of_interest, collapse = ",")), "Cell.Type", TRUE, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), ifelse(is.null(cell_types_of_interest), "NULL", paste(cell_types_of_interest, collapse = ",")), "Cell.Type", TRUE, "No error.")     
    }
    index <- index + 1 
  }
  for (feature_colname in feature_colname_options) {
    error_output <- tryCatch(
      expr = {
        calculate_cell_proportions3D(curr_spe,
                                     cell_types_of_interest = c("A", "B"),
                                     feature_colname = feature_colname,
                                     plot_image = TRUE)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", feature_colname, TRUE, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", feature_colname, TRUE, "No error.")     
    }
    index <- index + 1 
  }
      
  for (plot_image in plot_image_options) {
    error_output <- tryCatch(
      expr = {
        calculate_cell_proportions3D(curr_spe,
                                     cell_types_of_interest = c("A", "B"),
                                     feature_colname = "Cell.Type",
                                     plot_image = plot_image)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", "Cell.Type", plot_image, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", "Cell.Type", plot_image, "No error.")     
    }
    index <- index + 1 
  }
}

### Calculate entropy background ---------
error_df <- data.frame(matrix(nrow = 400, ncol = 4))
colnames(error_df) <- c("spe", "cell_type_of_interest", "feature_colname", "error_message")
index <- 1

cell_types_of_interest_options <- list(NULL, "A", c("A", "B"), c("A", "B", "C"), c("A", "B", "C", "D"), "Invalid")
feature_colname_options <- c("Cell.Type", "Invalid")

for (i in seq(length(spes))) {
  
  curr_spe <- spes[[i]]
  
  for (cell_types_of_interest in cell_types_of_interest_options) {
    error_output <- tryCatch(
      expr = {
        calculate_entropy_background3D(curr_spe,
                                       cell_types_of_interest = cell_types_of_interest,
                                       feature_colname = "Cell.Type")
      },
      error = function(e) {
        return(paste(e))
      }
    )
    error_df[index, ] <- c(paste("spe", i, sep = ""), ifelse(is.null(cell_types_of_interest), "NULL", paste(cell_types_of_interest, collapse = ",")), "Cell.Type", error_output)
    index <- index + 1 
  }
  for (feature_colname in feature_colname_options) {
    error_output <- tryCatch(
      expr = {
        calculate_entropy_background3D(curr_spe,
                                       cell_types_of_interest = c("A", "B"),
                                       feature_colname = feature_colname)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", feature_colname, error_output)
    index <- index + 1 
  }
}




### Calculate pairwise distances between cells -----------
pairwise_distances <- calculate_pairwise_distances_between_cell_types3D(chosen_spe,
                                                                        cell_types_of_interest = c("A", "B"),
                                                                        feature_colname = "Cell.Type",
                                                                        show_summary = TRUE,
                                                                        plot_image = TRUE)


### Calculate minimum distances between cells ------
error_df <- data.frame(matrix(nrow = 400, ncol = 6))
colnames(error_df) <- c("spe", "cell_type_of_interest", "feature_colname", "show_summary", "plot_image", "error_message")
index <- 1

cell_types_of_interest_options <- list(NULL, "A", c("A", "B"), c("A", "B", "C"), c("A", "B", "C", "D"), "Invalid")
feature_colname_options <- c("Cell.Type", "Invalid")
show_summary_options <- list(TRUE, FALSE, "Invalid")
plot_image_options <- list(TRUE, FALSE, "Invalid")

for (i in seq(length(spes))) {
  
  curr_spe <- spes[[i]]
  
  for (cell_types_of_interest in cell_types_of_interest_options) {
    error_output <- tryCatch(
      expr = {
        calculate_minimum_distances_between_cell_types3D(curr_spe,
                                                         cell_types_of_interest = cell_types_of_interest,
                                                         feature_colname = "Cell.Type",
                                                         show_summary = TRUE,
                                                         plot_image = TRUE)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), ifelse(is.null(cell_types_of_interest), "NULL", paste(cell_types_of_interest, collapse = ",")), "Cell.Type", TRUE, TRUE, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), ifelse(is.null(cell_types_of_interest), "NULL", paste(cell_types_of_interest, collapse = ",")), "Cell.Type", TRUE, TRUE, "No error.")     
    }
    index <- index + 1 
  }
  for (feature_colname in feature_colname_options) {
    error_output <- tryCatch(
      expr = {
        calculate_minimum_distances_between_cell_types3D(curr_spe,
                                                         cell_types_of_interest = c("A", "B"),
                                                         feature_colname = feature_colname,
                                                         show_summary = TRUE,
                                                         plot_image = TRUE)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", feature_colname, TRUE, TRUE, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", feature_colname, TRUE, TRUE, "No error.")     
    }
    index <- index + 1 
  }
  
  for (show_summary in show_summary_options) {
    error_output <- tryCatch(
      expr = {
        calculate_minimum_distances_between_cell_types3D(curr_spe,
                                                         cell_types_of_interest = c("A", "B"),
                                                         feature_colname = "Cell.Type",
                                                         show_summary = show_summary,
                                                         plot_image = TRUE)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", "Cell.Type", show_summary, TRUE, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", "Cell.Type", show_summary, TRUE, "No error.")     
    }
    index <- index + 1 
  }
  
  for (plot_image in plot_image_options) {
    error_output <- tryCatch(
      expr = {
        calculate_minimum_distances_between_cell_types3D(curr_spe,
                                                         cell_types_of_interest = c("A", "B"),
                                                         feature_colname = "Cell.Type",
                                                         show_summary = TRUE,
                                                         plot_image = plot_image)
      },
      error = function(e) {
        return(paste(e))
      }
    )
    
    if (is.character(error_output)) {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", "Cell.Type", TRUE, plot_image, error_output)
    } else {
      error_df[index, ] <- c(paste("spe", i, sep = ""), "A,B", "Cell.Type", TRUE, plot_image, "No error.")     
    }
    index <- index + 1 
  }
}




### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores3D(chosen_spe,
                                           reference_cell_types = c("A", "B"),
                                           target_cell_types = c("A", "B"),
                                           radius = 20,
                                           feature_colname = "Cell.Type",)
print(mixing_scores)

mixing_scores <- calculate_mixing_scores3D(chosen_spe,
                                           reference_cell_types = c("A", "B", "D"),
                                           target_cell_types = c("A", "B"),
                                           radius = 20,
                                           feature_colname = "Cell.Type",)
print(mixing_scores)

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(chosen_spe,
                                                             reference_cell_type = "A",
                                                             target_cell_type = "B",
                                                             radii = seq(10, 200, 10))

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(chosen_spe,
                                                             reference_cell_type = "B",
                                                             target_cell_type = "A",
                                                             radii = seq(10, 200, 10))


### Calculate cells in the neighbourhood
neighbourhood_cells <- calculate_cells_in_neighbourhood3D(chosen_spe,
                                                          reference_cell_type = "A",
                                                          target_cell_types = c("A", "B"),
                                                          radius = 30,
                                                          plot_image = F)

neighbourhood_cells_gradient <- calculate_cells_in_neighbourhood_gradient3D(chosen_spe,
                                                                            reference_cell_type = "A",
                                                                            target_cell_types = c("A", "B"),
                                                                            radii = seq(1, 30, 2),
                                                                            plot_image = T)


## Calculate cell proportions in the neighbourhood
neighbourhood_cell_proportions <- calculate_cells_in_neighbourhood_proportions3D(chosen_spe,
                                                                                 reference_cell_type = "A",
                                                                                 target_cell_types = c("A", "B"),
                                                                                 radius = 20)
print(neighbourhood_cell_proportions)

neighbourhood_cell_proportions_gradient <- calculate_cells_in_neighbourhood_proportions_gradient3D(chosen_spe,
                                                                                                   reference_cell_type = "A",
                                                                                                   target_cell_types = c("A", "B"),
                                                                                                   radii = seq(1, 50, 3))


### Calculate cross-K function
cross_K <- calculate_cross_K3D(chosen_spe,
                               reference_cell_type = "A",
                               target_cell_type = "B",
                               radius = 20)
print(cross_K)


cross_K_gradient <- calculate_cross_K_gradient3D(chosen_spe,
                                                 reference_cell_type = "A",
                                                 target_cell_type = "B",
                                                 radii = seq(1, 50))


### Calculate entropy
entropy_background <- calculate_entropy_background3D(chosen_spe,
                                                     cell_types_of_interest = c("A", "B"))

print(entropy_background)

entropy_result <- calculate_entropy3D(chosen_spe,
                                      radius = 20,
                                      reference_cell_type = "A",
                                      target_cell_types = c("A", "B"))





entropy_gradient <- calculate_entropy_gradient3D(chosen_spe,
                                                 reference_cell_type = "A",
                                                 target_cell_types = c("A", "B"),
                                                 radii = seq(1, 50, 2),
                                                 plot_image = TRUE)


### Using all_single_radius and all_gradient functions

all_single_radius_result <- calculate_all_single_radius_cc_metrics3D(chosen_spe, "A", c("A", "B"), 20)

all_gradient_result <- calculate_all_gradient_cc_metrics3D(chosen_spe, "A", c("A", "B"), seq(1, 50, 2))



### Calculate entropy grid metrics
entropy_grid_metrics <- calculate_entropy_grid_metrics3D(chosen_spe,
                                                         n_splits = 8,
                                                         cell_types_of_interest = c("A", "B", "B1"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(entropy_grid_metrics, "entropy")


### Calculate entropy prevalence
entropy_prevalence <- calculate_prevalence3D(entropy_grid_metrics,
                                             metric_colname = "entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

entropy_prevalence_gradient <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                               "entropy")

### Calculate spatial autocorrelation
entropy_spatial_autocorrelation <- calculate_spatial_autocorrelation3D(entropy_grid_metrics,
                                                                       metric_colname = "entropy",
                                                                       weight_method = "IDW")
print(entropy_spatial_autocorrelation)


### Calculate cell proportion grid metrics
cell_proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(chosen_spe,
                                                                         n_splits = 10,
                                                                         reference_cell_types = c("A"),
                                                                         target_cell_types = c("B"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(cell_proportion_grid_metrics, "proportion")


### Calculate cell proportion prevalence
cell_proportion_prevalence <- calculate_prevalence3D(cell_proportion_grid_metrics,
                                                     metric_colname = "proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)

cell_proportion_prevalence_gradient <- calculate_prevalence_gradient3D(cell_proportion_grid_metrics,
                                                                       metric_colname = "proportion")

## Calculate spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- calculate_spatial_autocorrelation3D(cell_proportion_grid_metrics,
                                                                               metric_colname = "proportion",
                                                                               weight_method = 0.10)
print(cell_proportion_spatial_autocorrelation)




spe_alpha_hull <- alpha_hull_clustering3D(chosen_spe, c("A", "B"), alpha = 15, minimum_cells_in_alpha_hull = 30)

plot_alpha_hull_clusters3D(spe_alpha_hull, c("A", "B", "C"), c("orange", "skyblue", "lightgray"))

alpha_hull_props <- calculate_cell_proportions_of_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster")

alpha_hull_min_distances <- calculate_minimum_distances_to_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster", 
                                                                      cell_types_inside_cluster = c("A", "B"),
                                                                      cell_types_outside_cluster = c("A"))

calculate_volume_of_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster")

calculate_center_of_clusters3D(spe_alpha_hull, "alpha_hull_cluster")

spe_alpha_hull <- calculate_border_of_clusters3D(spe_alpha_hull, 15, "alpha_hull_cluster")


spe_dbscan <- dbscan_clustering3D(chosen_spe, c("A", "B"), radius = 30, minimum_cells_in_radius = 20, minimum_cells_in_cluster = 30)

dbscan_props <- calculate_cell_proportions_of_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster")

dbscan_min_distances <- calculate_minimum_distances_to_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster", 
                                                                  cell_types_inside_cluster = c("A", "B"),
                                                                  cell_types_outside_cluster = c("A"))

calculate_volume_of_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster")

calculate_center_of_clusters3D(spe_dbscan, "dbscan_cluster")

spe_dbscan <- calculate_border_of_clusters3D(spe_dbscan, 6, "dbscan_cluster")


spe_grid <- grid_based_clustering3D(chosen_spe, cell_types_of_interest = c("A", "B"), n_splits = 10, minimum_cells_in_cluster = 30)

plot_grid_based_clusters3D(spe_grid, c("A", "B", "C"), c("orange", "skyblue", "lightgray"))

grid_props <- calculate_cell_proportions_of_clusters3D(spe_grid, cluster_colname = "grid_based_cluster")

grid_min_distances <- calculate_minimum_distances_to_clusters3D(spe_grid, cluster_colname = "grid_based_cluster", 
                                                                cell_types_inside_cluster = c("A", "B"),
                                                                cell_types_outside_cluster = c("A"))

calculate_volume_of_clusters3D(spe_grid, cluster_colname = "grid_based_cluster")

calculate_center_of_clusters3D(spe_grid, "grid_based_cluster")

spe_grid <- calculate_border_of_clusters3D(spe_grid, 8, "grid_based_cluster")



plot_cells3D(chosen_spe,
             plot_cell_types = c("A", "B", "C"),
             plot_colours = c("orange", "skyblue", "lightgray"))
