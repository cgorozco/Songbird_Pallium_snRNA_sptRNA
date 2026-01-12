library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr) # for stat_cor
library(EnhancedVolcano)
library(purrr)

plot_module_scores <- function(seuObj,
                               grouping = c("clusterName3", "BirdID"),
                               mod1, mod2, 
                               filterClusters) {
  
  # Fetch data
  allMods <- FetchData(seuObj, vars = c(grouping,
                                        paste0(mod1, "1"),
                                        paste0(mod2, "1")))
  
  # Aggregate averages
  allMods_avg <- allMods %>% 
    group_by(clusterName3, BirdID) %>% 
    mutate(
      Mod1Avg_byBird = mean(.data[[paste0(mod1, "1")]]),
      Mod2Avg_byBird = mean(.data[[paste0(mod2, "1")]])
    ) %>% 
    ungroup() %>% 
    group_by(clusterName3) %>% 
    mutate(
      Mod1Avg_byCluster = mean(Mod1Avg_byBird),
      Mod2ModAvg_byCluster = mean(Mod2Avg_byBird)
    ) %>% 
    ungroup() %>% 
    dplyr::select(clusterName3, Mod1Avg_byCluster, Mod2ModAvg_byCluster) %>% 
    filter(clusterName3 %in% filterClusters) %>% 
    unique()
  
  # Plot
  p <- ggplot(allMods_avg, aes(x = Mod1Avg_byCluster, y = Mod2ModAvg_byCluster, label = clusterName3)) +
    geom_point(color = "#010101", size = 3) +
    geom_smooth(method = "lm", se = T, 
               color = "darkgrey", linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(max.overlaps = 50) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
    theme_pubr() +
    #scale_x_continuous(limits = c(0, NA)) +
    #scale_y_continuous(limits = c(0, NA)) +
    labs(
      x = paste0(mod1, " Average Module Score"),
      y = paste0(mod2, " Average Module Score")
    ) + #coord_equal() +
    theme(
      axis.line = element_line(color = "#020201") 
    )

  return(p)
}


get_top_genes <- function(degs, 
                          direction = c("pos", "neg"), 
                          n = 500, 
                          outputFile = Tables_dir,
                          out_name = "top_genes") {
  direction <- match.arg(direction)
  
  if (direction == "pos") {
    genes <- degs %>% 
      dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% 
      arrange(desc(avg_log2FC)) %>% 
      slice_head(n = n) %>% 
      pull(gene)
  } else {
    genes <- degs %>% 
      dplyr::filter(avg_log2FC < 0 & p_val_adj < 0.05) %>% 
      arrange(avg_log2FC) %>% 
      slice_head(n = n) %>% 
      pull(gene)
  }
  
  # Save to CSV
  gene_table <- degs %>% filter(gene %in% genes)
  
  out_file <- paste0(outputFile, "/", out_name, ".csv")
  write.csv(gene_table, out_file, row.names = FALSE)
  
  return(genes)  
}


make_point_colors <- function(degs_to_plot,
                              pos_overlap_degs = pos_meso_degs,
                              neg_overlap_degs = neg_meso_degs,
                              color_TheseGenes = character(),
                              TheseGenes_color = "aquamarine3") {
  
  cols <- with(degs_to_plot, dplyr::case_when(
    # highlight custom genes if provided (if vector is length 0, condition is FALSE for all rows)
    rownames(degs_to_plot) %in% color_TheseGenes & p_val_adj < 0.05 ~ TheseGenes_color,
    
    # highlight overlaps
    rownames(degs_to_plot) %in% pos_overlap_degs & p_val_adj < 0.05 ~ "goldenrod3",
    rownames(degs_to_plot) %in% neg_overlap_degs & p_val_adj < 0.05 ~ "chartreuse4",
    
    # default DEG coloring
    avg_log2FC > 0 & p_val_adj < 0.05 ~ "brown3",
    avg_log2FC < 0 & p_val_adj < 0.05 ~ "royalblue",
    
    # everything else
    TRUE ~ "grey50"
  ))
  
  names(cols) <- rownames(degs_to_plot)
  cols
}

point_colors_simple <- function(degs_to_plot,
                              color_TheseGenes = character(),
                              TheseGenes_color = "aquamarine3") {
  
  cols <- with(degs_to_plot, dplyr::case_when(
    # highlight custom genes if provided (if vector is length 0, condition is FALSE for all rows)
    rownames(degs_to_plot) %in% color_TheseGenes & p_val_adj < 0.05 ~ TheseGenes_color,
    
    # default DEG coloring
    avg_log2FC > 0 & p_val_adj < 0.05 ~ "#F2727B",
    avg_log2FC < 0 & p_val_adj < 0.05 ~ "#6AA9CC",

  # everything else
    TRUE ~ "grey50"
  ))

  names(cols) <- rownames(degs_to_plot)
  cols
}


top_genes_label <- function(degs_to_plot, 
                            pos_fc = 15, pos_pVal = 20,
                            neg_fc = 15, neg_pVal = 20){
  
  degs_to_plot <- degs_to_plot %>% filter(!grepl("^LOC", gene))
  
  top_fc_pos <- degs_to_plot %>% 
    filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    arrange(desc(abs(avg_log2FC))) %>% 
    pull(gene) %>% unique() %>% 
    head(pos_fc)
  
  top_pVal_pos <- degs_to_plot %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(p_val_adj) %>% 
    pull(gene) %>% unique() %>% 
    head(pos_pVal)
  
  top_fc_neg <- degs_to_plot %>% 
    filter(avg_log2FC < 0 & p_val_adj < 0.05) %>% 
    arrange(desc(abs(avg_log2FC))) %>% 
    pull(gene) %>% unique() %>% 
    head(neg_fc)
  
  top_pVal_neg <- degs_to_plot %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(p_val_adj) %>% 
    pull(gene) %>% unique() %>% 
    head(neg_pVal)
  
  genes_to_label <- unique(c(top_fc_pos, top_pVal_pos,
                             top_fc_neg, top_pVal_neg))
  return(genes_to_label)
}

###===== code my into functions 




make_volcano <- function(degs_to_plot, out_file, title = "", asd_genes = character(), pos_meso_degs = character(), neg_meso_degs = character(), 
                         top_fc_n = 25, top_p_n = 35, labSize = 7, save = TRUE, width = 10, height = 12){
  point_colors <- make_point_colors(degs_to_plot, pos_meso_degs, neg_meso_degs, asd_genes)
  top_fc <- degs_to_plot %>% filter(p_val_adj < 0.005) %>% arrange(desc(abs(avg_log2FC))) %>% pull(gene) %>% unique() %>% head(top_fc_n)
  top_pVal <- degs_to_plot %>% arrange(p_val_adj) %>% pull(gene) %>% unique() %>% head(top_p_n)
  genes_to_label <- unique(c(top_fc, top_pVal))
  plt <- EnhancedVolcano(
    degs_to_plot,
    lab = degs_to_plot$gene,
    selectLab = genes_to_label,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    pCutoff = 0.05,
    FCcutoff = 0,
    cutoffLineCol = "#010101",
    drawConnectors = TRUE,
    arrowheads = FALSE,
    colCustom = point_colors,
    legendPosition = "none",
    gridlines.minor = FALSE,
    gridlines.major = FALSE,
    labSize = labSize,
    max.overlaps = 15
  ) + theme(axis.line = element_line(color = "#020201")) + ggtitle(title)
  if(save){
    ggsave(out_file, plot = plt, width = width, height = height)
  }
  plt
}

# Main analysis function (one call does what your repeated blocks did)
analyze_group <- function(seuObj,
                          ident1,                     # vector of clusters to treat as ident1 (e.g., nidoClusters or arcoClusters)
                          ident2,                     # vector of clusters to treat as ident2 (e.g., mesoClusters)
                          clustsToComp,               # song clusters to compare individually
                          name,                       # short name used in filenames and module names, e.g., "nido" or "arco"
                          meso_pattern = "Meso",      # pattern to identify mesopallium clusters
                          Figures_dir, Tables_dir,
                          asd_genes = character(),
                          assay_for_findMarkers = "SCTorigIdent3",
                          agg_assay = "RNA",
                          max_cells_per_ident = 1000,
                          top_n_genes = 500){
  # ensure character
  mesoClusters <- unique(seuObj$clusterName3[grepl(meso_pattern, seuObj$clusterName3)]) %>% as.character()
  # 1) DE between ident1 and ident2 to define module genes
  degs_between <- FindMarkers(seuObj, assay = assay_for_findMarkers,
                              ident.1 = ident1, ident.2 = ident2,
                              max.cells.per.ident = max_cells_per_ident,
                              recorrect_umi = FALSE)
  degs_between$gene <- rownames(degs_between)
  degs_between <- degs_between %>% filter(p_val_adj < 0.05 & !grepl("^LOC", gene))
  pos_genes <- get_top_genes(degs_between, "pos", n = top_n_genes)
  neg_genes <- get_top_genes(degs_between, "neg", n = top_n_genes)
  # dynamic module names so repeated runs don't clash
  mod_pos_name <- paste0(name, "_posMod")
  mod_neg_name <- paste0(name, "_negMod")
  seuObj <<- AddModuleScore(object = seuObj, features = list(neg_genes), name = mod_neg_name) # neg = "meso-like"
  seuObj <<- AddModuleScore(object = seuObj, features = list(pos_genes), name = mod_pos_name) # pos = ident1-like
  # plot module scores (uses your plot_module_scores function)
  mod_plot <- plot_module_scores(seuObj, mod1 = mod_neg_name, mod2 = mod_pos_name,
                                 filterClusters = c(mesoClusters, ident1, clustsToComp))
  ggsave(file.path(Figures_dir, paste0("moduleplot_", name, "SongGluts.svg")), mod_plot, width = 5, height = 5)
  
  # 2) pseudobulk (aggregate) for these clusters
  subset_clusters <- unique(c(mesoClusters, ident1, clustsToComp))
  seuObj_toCompare_degs <- subset(seuObj, subset = clusterName3 %in% subset_clusters)
  seuObj_agg <- AggregateExpression(seuObj_toCompare_degs,
                                    group.by = c("clusterName3", "BirdID"),
                                    return.seurat = TRUE)
  seuObj_agg$clusterName3 <- as.character(seuObj_agg$clusterName3)
  # make a dynamic category column 
  cat_col <- paste0(name, "_cat")
  seuObj_agg[[cat_col]] <- dplyr::case_when(
    seuObj_agg$clusterName3 %in% mesoClusters ~ "meso",
    seuObj_agg$clusterName3 %in% ident1 ~ name,
    seuObj_agg$clusterName3 %in% clustsToComp ~ "song",
    TRUE ~ seuObj_agg$clusterName3
  )
  # 3) mesos vs ident1 (pseudobulk DE)
  degs_meso <- FindMarkers(object = seuObj_agg, assay = agg_assay,
                           group.by = cat_col, ident.1 = "meso", ident.2 = name,
                           test.use = "DESeq2")
  degs_meso$gene <- rownames(degs_meso)
  write_csv(degs_meso[degs_meso$p_val_adj < 0.05,], file = file.path(Tables_dir, paste0("pseudoDegs_mesoV", name, ".csv")))
  pos_meso_degs <- degs_meso$gene[degs_meso$avg_log2FC > 0 & degs_meso$p_val_adj < 0.05 & !grepl("^LOC", degs_meso$gene)]
  neg_meso_degs <- degs_meso$gene[degs_meso$avg_log2FC < 0 & degs_meso$p_val_adj < 0.05 & !grepl("^LOC", degs_meso$gene)]
  
  # 4) song vs ident1 (pseudobulk DE) + volcano
  degs_song <- FindMarkers(object = seuObj_agg, assay = agg_assay, group.by = cat_col, ident.1 = "song", ident.2 = name, test.use = "DESeq2")
  degs_song$gene <- rownames(degs_song)
  write_csv(degs_song[degs_song$p_val_adj < 0.05,], file = file.path(Tables_dir, paste0("pseudoDegs_songGroupedVs", name, ".csv")))
  degs_to_plot <- degs_song %>% filter(!grepl("^LOC", gene))
  make_volcano(degs_to_plot,
               out_file = file.path(Figures_dir, paste0("volcano_Songvs", name, "_asdAqua.svg")),
               title = paste0("Song vs ", tools::toTitleCase(name), " (ASD in aqua)"),
               asd_genes = asd_genes,
               pos_meso_degs = pos_meso_degs,
               neg_meso_degs = neg_meso_degs,
               top_fc_n = 25, top_p_n = 35, labSize = 10, save = TRUE)
  
  # 5) per-cluster comparisons (song clusters individually vs ident1) and combined volcano grid
  degs_song_list <- list()
  for(cluster in clustsToComp){
    degs_cluster <- FindMarkers(object = seuObj_agg, assay = agg_assay,
                                group.by = paste0(name, "_cat"),
                                ident.1 = cluster, ident.2 = name,
                                test.use = "DESeq2")
    degs_cluster$gene <- rownames(degs_cluster)
    write_csv(degs_cluster[degs_cluster$p_val_adj < 0.05,], file = file.path(Tables_dir, paste0("pseudoDegs_", name, "Vs", cluster, ".csv")))
    degs_song_list[[cluster]] <- degs_cluster
  }
  
  # build volcano list (one plot per song cluster)
  plot_list <- lapply(clustsToComp, function(cluster){
    degs_to_plot <- degs_song_list[[cluster]] %>% filter(!grepl("^LOC", gene))
    point_colors <- make_point_colors(degs_to_plot, pos_meso_degs, neg_meso_degs, asd_genes)
    top_fc_pos <- degs_to_plot %>% filter(avg_log2FC > 1 & p_val_adj < 0.005) %>% arrange(desc(avg_log2FC)) %>% pull(gene) %>% unique() %>% head(7)
    top_pVal_pos <- degs_to_plot %>% arrange(p_val_adj) %>% pull(gene) %>% unique() %>% head(15)
    genes_to_label <- unique(c(top_fc_pos, top_pVal_pos))
    EnhancedVolcano(
      degs_to_plot,
      lab = degs_to_plot$gene,
      selectLab = genes_to_label,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      pCutoff = 0.05,
      FCcutoff = 0,
      cutoffLineCol = "#010101",
      drawConnectors = TRUE,
      arrowheads = FALSE,
      colCustom = point_colors,
      legendPosition = "none",
      gridlines.minor = FALSE,
      gridlines.major = FALSE,
      labSize = 7,
      max.overlaps = 15
    ) + theme(axis.line = element_line(color = "#020201")) + ggtitle(cluster)
  })
  
  combined_plot <- wrap_plots(plot_list, ncol = min(4, length(plot_list)))
  ggsave(file.path(Figures_dir, paste0("volcanosCombined_", name, "SongGluts_upNdownMeso.svg")), combined_plot, width = 12, height = 7)
  
  # return a list of important objects for further inspection
  list(
    degs_between = degs_between,
    pos_genes = pos_genes,
    neg_genes = neg_genes,
    degs_meso = degs_meso,
    degs_song = degs_song,
    degs_song_list = degs_song_list,
    module_plot = mod_plot,
    combined_volcano = combined_plot
  )
}



