library(patchwork)

# CREATE UMAP WITH ONLY CLUSTER HIGHLIGHTED AT A TIME AND SAVE INTO SINGLE FILE
highlight_subclusters <- function(seuObj, clustering, filter = F, pattern, reduction = "umap.integrated.rpca",
                                  Figures_dir = Figures_dir, ncol = 5, width = 30, height = 45) {
  
  IDs <- unique(seuObj[[clustering]][[1]])
  
  if (filter){
    IDs <- IDs[grepl(pattern, IDs)]
  }

  Idents(seuObj) <- clustering
  plot_list <- list()
  for (id in IDs) {
    p <- DimPlot(seuObj, reduction = reduction,
                 cells.highlight = CellsByIdentities(seuObj, idents = id), order = id) + coord_equal()
    plot_list[[as.character(id)]] <- p
  }
  combined_plot <- wrap_plots(plot_list, ncol = ncol)
  ggsave(filename = paste0(Figures_dir, "/umap_HiglightBy", clustering, ".tiff"),  plot = combined_plot,
         bg = "white", width = width, height = height)
  
  # Return the combined plot
  return(combined_plot)
  
}


# CREATE HEATMAP OF AVERAGE LOG EXPRESSION OF SELECTED GENES
#   AverageExpression() MUST HAVE ALREADY BEEN RUN
avgExpr_heatmap <- function(clustAvgUse, genesPlot,
                            clustersUse, clusterOrder = NULL,
                            colOpt = "magma", scaleOpt = F, wht2blk = F) {
  
  clustAvg_plot <- FetchData(object = clustAvgUse, vars = genesPlot, layer = "data")
  
  df_long <- clustAvg_plot %>%
    as.data.frame() %>%
    rownames_to_column(var = clustersUse) %>%
    pivot_longer(cols = -!!sym(clustersUse), names_to = "Gene", values_to = "Avg")
  
  if(!is.null(clusterOrder)){
    df_long[[clustersUse]] <- factor(df_long[[clustersUse]], levels = gsub("_", "-", clusterOrder))
  }
  
  df_long$Gene <- factor(df_long$Gene, levels = genesPlot)

  CanMrks_heat <- ggplot(df_long, aes(x = Gene, y = !!sym(clustersUse), fill = Avg)) +
    geom_tile() +
    scale_fill_viridis_c(option = colOpt) +
    theme_minimal() +
    labs(fill = "Avg. Log Expr.") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    coord_equal()
  
  # Plot white to black expression
  if (wht2blk == T){
    CanMrks_heat <- ggplot(df_long, aes(x = Gene, y = !!sym(clustersUse), fill = Avg)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      labs(fill = "Avg. Log Expr.") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      coord_equal()
  }


  # Plot with scaled values
  if (scaleOpt == T){
    df_long_scaled <- df_long %>%
      group_by(Gene) %>%
      mutate(Avg_scaled = scale(Avg)) %>%
      ungroup()
    
    CanMrks_heat <- ggplot(df_long_scaled, aes(x = Gene, y = !!sym(clustersUse), fill = Avg_scaled)) +
      geom_tile() +
      scale_fill_viridis_c(option = colOpt) +
      theme_minimal() +
      labs(fill = "Scaled Avg. \nLog Expr.") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      coord_equal()
  }

  return(CanMrks_heat)
  
}




# CREATE HEATMAP FOR AVERAGE PREDICTION SCORES
predicScore_heatmap <- function(zf_meta, predic_meta, filterZF, filterZF_col, groupUse, otherFilter,
                                zf_order = NULL, other_order = NULL){
  
  metaComb <- full_join(
    rownames_to_column(zf_meta, var = "rowname"),
    rownames_to_column(predic_meta, var = "rowname"),
    by = "rowname"
  )

  heatToPlot <- metaComb %>%
    filter(grepl(filterZF, !!sym(filterZF_col))) %>%
    select(c(!!sym(filterZF_col), !!sym(groupUse)), contains(otherFilter)) %>% 
    group_by(!!sym(groupUse)) %>%
    summarize(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
    column_to_rownames(groupUse) %>%
    select(where(~ any(. > 0.05)))
  
  colnames(heatToPlot) <- gsub("prediction.score.", "", colnames(heatToPlot))

  heatToPlot_long <- as.data.frame(heatToPlot) %>%
    rownames_to_column(var = "Cluster") %>%
    pivot_longer(cols = -Cluster, names_to = "PredicClass", values_to = "Score")
  
  if (!is.null(zf_order)){
    heatToPlot_long$Cluster <- factor(heatToPlot_long$Cluster, levels = zf_order)
  }
  
  if (!is.null(other_order)){
    heatToPlot_long$PredicClass <- factor(heatToPlot_long$PredicClass, levels = other_order)
  }
  
  heatPlot_gg <- ggplot(heatToPlot_long, aes(x = PredicClass, y = Cluster, fill = Score)) +
    geom_tile(color = "lightgray", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    theme_minimal() +
    coord_equal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

  return(heatPlot_gg)
  # 
  # heatPlot <- pheatmap(heatToPlot, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("white", "firebrick"))(100))
  # return(heatPlot)
  
}











































