make_bubble_plot <- function(cellchat, df.net, 
                             clustsToComp,        # vector of comparison clusters
                             cluster_pattern,     # regex for background clusters
                             targets.use,         # vector of target clusters
                             filename,            # filename to save plot
                             Figures_dir) {
  
  # if cluster_pattern is a vector, use %in%, otherwise treat it as regex
  if (length(cluster_pattern) > 1) {
    clusters <- cellchat@meta[["clusterName3"]][cellchat@meta[["clusterName3"]] %in% cluster_pattern] %>%
      unique() %>% as.character()
  } else {
    clusters <- cellchat@meta[["clusterName3"]][grepl(cluster_pattern, cellchat@meta[["clusterName3"]])] %>%
      unique() %>% as.character()
  }
  
  # filter informative pathways
  df_net_use <- df.net %>% 
    filter(source %in% c(clusters, clustsToComp) & target %in% targets.use) %>% 
    group_by(interaction_name) %>% 
    filter(any(source %in% clustsToComp) & !any(source %in% clusters)) %>% 
    ungroup()
  
  # collect pairs
  pair_df <- df_net_use %>% 
    distinct(interaction_name, ligand) %>%             # keep ligand info
    arrange(factor(ligand, levels = unique(df_net_use$ligand))) %>%  # order by ligand
    mutate(interaction_name = factor(interaction_name, levels = interaction_name)) %>% 
    select(interaction_name) %>%                      # keep only interaction_name
    as.data.frame()
  
  # make plot
  p <- netVisual_bubble(cellchat, 
                        sources.use = c(clustsToComp, clusters),
                        targets.use = targets.use,
                        pairLR.use = pair_df,
                        color.heatmap = "viridis",
                        remove.isolate = FALSE) + 
    coord_equal()

  # save plot
  ggsave(filename = file.path(Figures_dir, filename),
         plot = p,
         width = 4, height = 7)
  
  # return both
  return(list(pair_df = pair_df, plot = p))
}
