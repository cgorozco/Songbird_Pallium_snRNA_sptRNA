source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/correlations.R")

library(reshape2)
library(ComplexHeatmap)

generate_SpearmanCorrHeatmap_wSig <- function(exprAvg_filt_zf, exprAvg_filt_other,
                                              renameClusters = F, 
                                              glut.rename, zf.levels, 
                                              other.levels = NULL,
                                              nrep = 500, 
                                              dropSelf = F,
                                              minSetTo,
                                              FigName, Figures_dir, widthUse = 6) {
  # Log-transform and adjust matrices
  mat_a <- log1p(exprAvg_filt_zf) + 0.1
  mat_b <- log1p(exprAvg_filt_other) + 0.1

  # Calculate correlations
  obj_cor <- specificity_correlate(mat_a, mat_b, method = "spearman")
  
  if (dropSelf) {
    diag(obj_cor) <- NA
  }
  
  # Recode rownames using `glut.rename` and order them based on `zf.levels`
  if (renameClusters == T){
    rownames(obj_cor) <- recode(rownames(obj_cor), !!!glut.rename)
    obj_cor <- obj_cor[order(factor(rownames(obj_cor), levels = zf.levels)), ]
  }
  
  # Or just reorder 
  if (!is.null(zf.levels)){
    obj_cor <- obj_cor[order(factor(rownames(obj_cor), levels = rev(zf.levels))), ]
  }
  
  if (!is.null(other.levels)){
    obj_cor <- obj_cor[, order(factor(colnames(obj_cor), levels = other.levels))]
  }
  
  # Shuffle correlations and compute statistics
  obj_cor_shuf <- map(1:nrep, function(i) {
    mat_a_cur <- mat_a
    mat_b_cur <- mat_b
    
    # Shuffle the rownames of mat_a_cur
    rownames(mat_a_cur) <- sample(rownames(mat_a_cur))
    
    # Perform the correlation
    obj_cor <- specificity_correlate(mat_a_cur, mat_b_cur)
    
    # Reshape the correlation result into a dataframe
    obj_cor_df <- melt(obj_cor)
    colnames(obj_cor_df) <- c("celltype", "compare", "value")
    obj_cor_df <- obj_cor_df %>% filter(!is.na(value))
  
    return(obj_cor_df)
    }) %>% bind_rows()

  # Summarize shuffled results
  obj_cor_shuf_stat <- obj_cor_shuf %>%
    group_by(celltype, compare) %>%
    summarize(value_mean = mean(value),
              value_sd = sd(value),
              value_q99 = quantile(value, .99),
              value_q95 = quantile(value, .95),
              value_q995 = quantile(value, .995))
  
  # Melt the original correlation object for annotation
  obj_cor_df <- melt(obj_cor)
  colnames(obj_cor_df) <- c("celltype", "compare", "value")
  
  # Recode celltype in shuffled stats and join with original correlations
  if (renameClusters == T){
    obj_cor_shuf_stat$celltype <- factor(recode(obj_cor_shuf_stat$celltype, !!!glut.rename), levels = zf.levels)
  }
  
  obj_cor_df <- obj_cor_df %>%
    left_join(obj_cor_shuf_stat, by = c("celltype", "compare"))
  
  # Add significance label
  obj_cor_df <- obj_cor_df %>%
    mutate(sig_label = if_else(value > value_q995, "*", ""))
  obj_cor_sig_adj_label <- acast(obj_cor_df, celltype ~ compare, value.var = "sig_label")
  
  # Define the star annotation function
  star_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(obj_cor_sig_adj_label[i, j]) && obj_cor_sig_adj_label[i, j] != "") {
      grid.text(obj_cor_sig_adj_label[i, j], x, y, gp = gpar(fontsize = widthUse*1.5))
    }
  }
  
  # set color
  col_fun <- colorRampPalette(c("white", "#aa7c00"))(100)
  
  # set min to 0
  if (!is.null(minSetTo)){
    obj_cor[obj_cor < 0] <- minSetTo
  }


  # Create heatmap with clustering
  hm <- Heatmap(obj_cor, 
                cell_fun = star_fun, 
                cluster_rows = TRUE, 
                cluster_columns = TRUE, 
                show_row_names = TRUE,  
                show_column_names = TRUE,
                col = col_fun,
                width = unit(ncol(obj_cor) * 0.5, "cm"), 
                height = unit(nrow(obj_cor) * 0.5, "cm")
                )
  
  # Save clustered heatmap as TIFF
  heightUse <- widthUse*0.8
  tiff(file.path(Figures_dir, paste0("/heatmap_", FigName, ".tiff")), width = widthUse, height = heightUse, units = "in", res = 300)
  draw(hm)
  dev.off()
  
  # Create heatmap without clustering
  hm <- Heatmap(obj_cor, 
                cell_fun = star_fun, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                show_row_names = TRUE,  
                show_column_names = TRUE,
                col = col_fun,
                width = unit(ncol(obj_cor) * 0.5, "cm"), 
                height = unit(nrow(obj_cor) * 0.5, "cm"),
                rect_gp = gpar(col = "gray", lwd = 1)
  )
  
  return(list(
    hm = hm,
    sig_labels = obj_cor_sig_adj_label,
    cor_matrix = obj_cor
  ))
  
  }












###===== with a floor of 0
generate_SpearmanCorrHeatmap_wFloorNSig <- function(exprAvg_filt_zf, exprAvg_filt_other,
                                              renameClusters = TRUE, glut.rename, zf.levels, other.levels,
                                              nrep = 100, 
                                              FigName, Figures_dir) {
  
  # Log-transform and adjust matrices
  mat_a <- log1p(exprAvg_filt_zf) + 0.1
  mat_b <- log1p(exprAvg_filt_other) + 0.1
  
  # Calculate correlations
  obj_cor <- specificity_correlate(mat_a, mat_b, method = "spearman")
  
  # Set values lower than 0 to 0
  obj_cor <- pmax(obj_cor, 0)
  
  # Recode rownames using `glut.rename` and order them based on `zf.levels`
  if (renameClusters) {
    rownames(obj_cor) <- recode(rownames(obj_cor), !!!glut.rename)
    obj_cor <- obj_cor[order(factor(rownames(obj_cor), levels = zf.levels)), ]
    obj_cor <- obj_cor[, order(factor(colnames(obj_cor), levels = other.levels))]
  }
  
  # Shuffle correlations and compute statistics
  obj_cor_shuf <- map(1:nrep, function(i) {
    mat_a_cur <- mat_a
    mat_b_cur <- mat_b
    
    # Shuffle the rownames of mat_a_cur
    rownames(mat_a_cur) <- sample(rownames(mat_a_cur))
    
    # Perform the correlation
    obj_cor <- specificity_correlate(mat_a_cur, mat_b_cur)
    
    # Set values lower than 0 to 0
    obj_cor <- pmax(obj_cor, 0)
    
    # Reshape the correlation result into a dataframe
    obj_cor_df <- melt(obj_cor)
    colnames(obj_cor_df) <- c("celltype", "compare", "value")
    
    return(obj_cor_df)
  }) %>% bind_rows()
  
  # Summarize shuffled results
  obj_cor_shuf_stat <- obj_cor_shuf %>%
    group_by(celltype, compare) %>%
    summarize(value_mean = mean(value),
              value_sd = sd(value),
              value_q99 = quantile(value, .99),
              value_q95 = quantile(value, .95))
  
  # Melt the original correlation object for annotation
  obj_cor_df <- melt(obj_cor)
  colnames(obj_cor_df) <- c("celltype", "compare", "value")
  
  # Recode celltype in shuffled stats and join with original correlations
  if (renameClusters) {
    obj_cor_shuf_stat$celltype <- factor(recode(obj_cor_shuf_stat$celltype, !!!glut.rename), levels = zf.levels)
  }
  obj_cor_df <- obj_cor_df %>%
    left_join(obj_cor_shuf_stat, by = c("celltype", "compare"))
  
  # Add significance label
  obj_cor_df <- obj_cor_df %>%
    mutate(sig_label = if_else(value > value_q95, "*", ""))
  obj_cor_sig_adj_label <- acast(obj_cor_df, celltype ~ compare, value.var = "sig_label")
  
  # Define the star annotation function
  star_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(obj_cor_sig_adj_label[i, j]) && obj_cor_sig_adj_label[i, j] != "") {
      grid.text(obj_cor_sig_adj_label[i, j], x, y, gp = gpar(fontsize = 10))
    }
  }
  
  col_fun <- colorRampPalette(c("white", "red"))(100)

  # Create heatmap with clustering
  hm <- Heatmap(obj_cor, 
                cell_fun = star_fun, 
                cluster_rows = TRUE, 
                cluster_columns = TRUE, 
                show_row_names = TRUE,  
                show_column_names = TRUE,
                col = col_fun)  # Apply color scale
  
  # Save clustered heatmap as TIFF
  tiff(file.path(Figures_dir, paste0("/heatmap_", FigName, ".tiff")), width = 10, height = 8, units = "in", res = 300)
  draw(hm)
  dev.off()
  
  # Create heatmap without clustering
  hm <- Heatmap(obj_cor, 
                cell_fun = star_fun, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                show_row_names = TRUE,  
                show_column_names = TRUE,
                col = col_fun)  # Apply color scale
  
  # Save non-clustered heatmap as TIFF
  tiff(file.path(Figures_dir, paste0("/heatmap_", FigName, "_noDendo.tiff")), width = 10, height = 8, units = "in", res = 300)
  draw(hm)
  dev.off()
}

# Example call to the function:
# generate_SpearmanCorrHeatmap_wSig(exprAvg_filt_zf, exprAvg_filt_other, 
#                                   glut.rename = custom_glut_rename, 
#                                   zf.levels = custom_glut_levels, 
#                                   other.levels = custom_other_levels, 
#                                   nrep = 100, 
#                                   FigName = "my_figure", 
#                                   Figures_dir = "/path/to/figures")