# Make in situ searching probe
create_inSituTest_plots <- function(degs_df, seu_obj, cluster_pattern, logFC_thresh, pct2_thresh, include_LOC, num_genes, fig_dir, fig_name) {
  # Filter based on whether to include LOC genes or not
  if (include_LOC) {
    probeTest <- degs_df %>%
      filter(grepl(cluster_pattern, cluster) & avg_log2FC > logFC_thresh & pct.2 < pct2_thresh) %>%
      arrange(desc(avg_log2FC)) %>%
      head(num_genes)
    loc_label <- "_wLoc"
  } else {
    probeTest <- degs_df %>%
      filter(grepl(cluster_pattern, cluster) & avg_log2FC > logFC_thresh & pct.2 < pct2_thresh & !grepl("^LOC", gene)) %>%
      arrange(desc(avg_log2FC)) %>%
      head(num_genes)
    loc_label <- ""
  }
  
  # Create and save DotPlot
  DotPlot(seu_obj, features = c(probeTest$gene)) + RotatedAxis()
  ggsave(filename = paste0(fig_dir, "/", fig_name, "_dotplot", loc_label, ".tiff"), bg = "white", width = 6, height = 10)
  
  # Create and save FeaturePlot
  FeaturePlot(seu_obj, features = c(probeTest$gene), coord.fixed = TRUE)
  ggsave(filename = paste0(fig_dir, "/", fig_name, "_featplot", loc_label, ".tiff"), bg = "white", width = 15, height = 12)
}

# # Example usage
# create_inSituTest_plots(allCells_degs, seuObj, "_11", 1, 0.1, include_LOC = F, num_genes = 9, fig_dir = Figures_dir, fig_name = "glut11")
# create_inSituTest_plots(allCells_degs, seuObj, "_11", 1, 0.1, include_LOC = T, num_genes = 9, fig_dir = Figures_dir, fig_name = "glut11")
