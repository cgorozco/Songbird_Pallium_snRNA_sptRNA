SmoothSpatialIdents4Merge <- function(seurat_obj, radius = NULL) {
  # Store original identities
  seurat_obj$original_idents <- Idents(seurat_obj)
  
  # Get unique slices identifiers
  slices_ids <- unique(seurat_obj$orig.ident)
  
  # Create list to store processed slices
  processed_slices <- list()
  
  # Process each slice
  for (slice_id in slices_ids) {
    # Subset the seurat object for this slice
    slice_subset <- subset(seurat_obj, subset = orig.ident == slice_id)
    
    # Apply smoothing to this slice
    slice_subset <- SmoothSpatialIdents(slice_subset, radius = radius)
    
    # Store processed slice
    processed_slices[[as.character(slice_id)]] <- slice_subset
  }
  
  # Merge processed slices back together
  # Start with the first slice
  merged_obj <- processed_slices[[1]]
  
  # Add remaining slices if they exist
  if (length(processed_slices) > 1) {
    for (i in 2:length(processed_slices)) {
      merged_obj <- merge(merged_obj, processed_slices[[i]])
    }
  }
  
  return(merged_obj)
}



SmoothSpatialIdents <- function(seurat_obj, radius = NULL) {
  # Get spatial coordinates
  coords <- GetTissueCoordinates(seurat_obj)
  
  # If radius is not specified, calculate based on average nearest neighbor distance
  if (is.null(radius)) {
    distances <- dist(coords)
    radius <- mean(distances) * 0.1
  }
  
  # Get current identities
  current_idents <- Idents(seurat_obj)
  new_idents <- current_idents
  
  # For each spot
  for (i in 1:nrow(coords)) {
    # Calculate distances to all other spots
    distances <- sqrt((coords[i,1] - coords[,1])^2 + (coords[i,2] - coords[,2])^2)
    
    # Find neighbors within radius
    neighbors <- which(distances <= radius & distances > 0)
    
    if (length(neighbors) > 0) {
      # Get identities of neighbors
      neighbor_idents <- current_idents[neighbors]
      
      # Count frequency of each identity in neighborhood
      ident_counts <- table(neighbor_idents)
      
      # If current cell's identity is different from the majority
      current_cell_ident <- current_idents[i]
      most_common_ident <- names(which.max(ident_counts))
      
      if (current_cell_ident != most_common_ident) {
        # Check if current identity is minority in neighborhood
        current_ident_count <- sum(neighbor_idents == current_cell_ident)
        most_common_count <- max(ident_counts)
        
        # Reassign if current identity is minority
        if (current_ident_count < most_common_count) {
          new_idents[i] <- most_common_ident
        }
      }
    }
  }
  
  # Assign new identities
  Idents(seurat_obj) <- new_idents
  
  # Return modified object
  return(seurat_obj)
}
