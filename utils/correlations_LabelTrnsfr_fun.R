run_label_transfer_shuffle_heatmap <- function(
    heatPlotUse = heatToPlot,
    gg_seuObj,
    gg_shuffleThis,
    zf_seuObj,
    anchors = anchors_GGtoZF,
    nrep = 300,
    seed = 42,
    col_use = colorRampPalette(c("white", "#aa7c00"))(100),
    RDSfiles,
    rds_name,
    reuse_shuffles = TRUE,
    force_rerun = FALSE
) {
  
  rds_path <- file.path(RDSfiles, rds_name)
  
  ## -------------------------
  ## SHUFFLE (OR LOAD)
  ## -------------------------
  if (reuse_shuffles && file.exists(rds_path) && !force_rerun) {
  
    message("Loading existing shuffled label-transfer results:")
    message("  ", rds_path)
    
    obj_cor_shuf <- readRDS(rds_path)
    
  } else {
    
    message("Running label-transfer shuffles (n = ", nrep, ")")
    
    set.seed(seed)
    
    obj_cor_shuf <- map(1:nrep, function(i) {
      
      Gg_tmp   <- gg_seuObj
      glut_tmp <- zf_seuObj
    
      Shuff_data <- FetchData(Gg_tmp, vars = gg_shuffleThis)
      Shuff_data$Shuffle1 <- sample(Shuff_data[[gg_shuffleThis]])
  
      Gg_tmp <- AddMetaData(Gg_tmp, Shuff_data$Shuffle1, col.name = "Shuffle1")
      
      predictions_shuff <- TransferData(
        anchorset = anchors_GGtoZF,
        refdata = Gg_tmp$Shuffle1,
        weight.reduction = "cca"
      )
  
      glut_tmp <- AddMetaData(glut_tmp, metadata = predictions_shuff)
      
      zf_meta <- FetchData(glut_tmp, vars = "clusterName3")
      
      metaComb <- full_join(
        rownames_to_column(zf_meta, var = "rowname"),
        rownames_to_column(predictions_shuff, var = "rowname"),
        by = "rowname"
      ) %>%
        select(-rowname, -predicted.id, -prediction.score.max)
      
      colnames(metaComb) <- gsub("prediction.score.", "", colnames(metaComb))
      
      metaComb_avg <- metaComb %>%
        group_by(clusterName3) %>%
        summarize(across(where(is.numeric), mean, na.rm = TRUE))
      
      obj_cor_df <- melt(metaComb_avg)
      colnames(obj_cor_df) <- c("celltype", "compare", "value")
      
      obj_cor_df %>% filter(!is.na(value))
    }) %>% bind_rows()
    
    saveRDS(obj_cor_shuf, rds_path)
  }
  
  ## -------------------------
  ## SHUFFLE STATS
  ## -------------------------
  obj_cor_shuf_stat <- obj_cor_shuf %>%
    group_by(celltype, compare) %>%
    summarize(
      value_mean = mean(value),
      value_sd   = sd(value),
      value_q99  = quantile(value, .99),
      value_q95  = quantile(value, .95),
      value_q995 = quantile(value, .995),
      value_q999 = quantile(value, .999),
      value_q9995 = quantile(value, .9995),
      .groups = "drop"
    )
  
  ## -------------------------
  ## REAL DATA + MATRICES
  ## -------------------------
  obj_cor <- rownames_to_column(heatPlotUse, var = "celltype")
  
  obj_cor_df <- melt(obj_cor)
  colnames(obj_cor_df) <- c("celltype", "compare", "value")
  
  obj_cor_df <- obj_cor_df %>%
    left_join(obj_cor_shuf_stat, by = c("celltype", "compare")) %>%
    mutate(sig_label = if_else(value > value_q9995, "*", ""))
  
  cor_mat <- obj_cor %>%
    column_to_rownames("celltype") %>%
    as.matrix()
  
  sig_label_mat <- acast(
    obj_cor_df,
    celltype ~ compare,
    value.var = "sig_label"
  )
  
  sig_label_mat <- sig_label_mat[
    rownames(cor_mat),
    colnames(cor_mat)
  ]
  
  ## -------------------------
  ## HEATMAP
  ## -------------------------
  star_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(sig_label_mat[i, j]) && sig_label_mat[i, j] != "") {
      grid.text("*", x, y, gp = gpar(fontsize = 9))
    }
  }
  
  col_fun <- col_use
  
  hm <- Heatmap(
    cor_mat,
    cell_fun = star_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    col = col_fun,
    width = unit(ncol(cor_mat) * 0.5, "cm"),
    height = unit(nrow(cor_mat) * 0.5, "cm"),
    rect_gp = gpar(col = "gray", lwd = 1)
  )
  
  plot_out <- as.ggplot(hm) + coord_equal()
  
  ## -------------------------
  ## RETURN
  ## -------------------------
  list(
    plot = plot_out,
    cor_mat = cor_mat,
    sig_label_mat = sig_label_mat,
    cor_df = obj_cor_df,
    shuffles = obj_cor_shuf
  )
}

