# SET UP
library(Seurat)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(ggrastr)
library(viridis)
library(scCustomize)

source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/correlations.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/spearmanCorr_heatmap.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/correlations_LabelTrnsfr_fun.R")

setwd("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/13_chicken")
root_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/13_chicken"

RDSfiles <- paste0(root_dir, "/RDS_files")
dir.create(RDSfiles, showWarnings = F)

Figures_dir <- paste0(root_dir, "/Figures/gluts")
dir.create(Figures_dir, showWarnings = F)

Tables_dir <- paste0(root_dir, "/Tables")
dir.create(Tables_dir, showWarnings = F)

# LOAD DATA
seuObj <- readRDS("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")
Gg_adult <- readRDS(file = paste0(RDSfiles, "/ChkZaremba25_wOrthos_sct.rds"))

Gg_adult_exc <- subset(Gg_adult, subset = tg_class == "Excitatory neurons" & anno_level_3 != "Ex_Pre")
Idents(Gg_adult_exc) <- "anno_level_3"

seuObj_glut <- subset(seuObj, subset = cn3_class %in% "Glutamatergic" & cn3_subclass != "Neurogenic")

re_do <- F
if (re_do == T) {
  Gg_adult_exc <- SCTransform(Gg_adult_exc,
    assay = "RNA_ortho", new.assay.name = "SCTexc",
    variable.features.n = 6000
  )
  saveRDS(Gg_adult_exc, file = paste0(RDSfiles, "/gg_Exc_SCT_20251020.rds"))

  seuObj_glut <- SCTransform(seuObj_glut,
    assay = "RNA", new.assay.name = "SCTexc",
    variable.features.n = 6000
  )
  saveRDS(seuObj_glut, file = paste0(RDSfiles, "/zf_Exc_SCT_20251020.rds"))
}

Gg_adult_exc <- readRDS(file = paste0(RDSfiles, "/gg_Exc_SCT_20251020.rds"))
seuObj_glut <- readRDS(file = paste0(RDSfiles, "/zf_Exc_SCT_20251020.rds"))

varFeats_AllCells <- intersect(Gg_adult@assays[["SCT_ortho"]]@var.features, seuObj@assays[["SCTorigIdent3"]]@var.features)
varFeats_GlutCells <- intersect(Gg_adult_exc@assays[["SCTexc"]]@var.features, seuObj_glut@assays[["SCTexc"]]@var.features)
varFeats_both <- union(varFeats_AllCells, varFeats_GlutCells)

redo2 <- F
if (redo2) {
  anchors_GGtoZF <- FindTransferAnchors(
    reference = Gg_adult_exc, query = seuObj_glut,
    features = varFeats_both,
    reference.assay = "SCTexc",
    query.assay = "SCTexc",
    normalization.method = "SCT",
    recompute.residuals = F,
    reduction = "cca"
  )

  # saveRDS(anchors_GGtoZF, file = paste0(RDSfiles, "/anchors_Exc_rpca_SCTorthoToOrigIdent3_20251021.rds"))
}

anchors_GGtoZF <- readRDS(file = paste0(RDSfiles, "/anchors_Exc_cca_SCTorthoToOrigIdent3_20251021.rds"))


predictions_join_ann3 <- TransferData(
  anchorset = anchors_GGtoZF, refdata = Gg_adult_exc$anno_level_3,
  weight.reduction = "cca"
)
seuObj_glut <- AddMetaData(seuObj_glut, metadata = predictions_join_ann3)

# VISUALIZE
zf_meta <- FetchData(seuObj_glut, vars = c("clusterName3"))

metaComb <- full_join(
  rownames_to_column(zf_meta, var = "rowname"),
  rownames_to_column(predictions_join_ann3, var = "rowname"),
  by = "rowname"
)

heatToPlot <- metaComb %>%
  filter(grepl("Glut", clusterName3)) %>%
  select(clusterName3, contains("Ex")) %>%
  # select(-prediction.score.Ex_Pre) %>%
  group_by(clusterName3) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames("clusterName3") %>%
  select(where(~ any(. > 0.00)))

colnames(heatToPlot) <- gsub("prediction.score.", "", colnames(heatToPlot))

heatToPlot <- heatToPlot[order(factor(rownames(heatToPlot), levels = rev(glutLvls_byRgnPct))), ]
heatToPlot <- heatToPlot[, order(factor(colnames(heatToPlot), levels = gg_levels_lvl3_v2))] # gg_levels_lvl3_v2

glut_corr_heatmap <- pheatmap(heatToPlot,
  cellwidth = 15, cellheight = 15,
  cluster_rows = F, cluster_cols = F,
  color = colorRampPalette(c("white", "#aa7c00", "#6b3e08"))(150)
) # colorRampPalette(c("white", "#aa7c00"))(100)) / viridis(100, option = "viridis") / rev(viridis(100, option = "mako"))

glut_corr_heatmap
p1 <- as.ggplot(glut_corr_heatmap) + coord_equal()

res <- run_label_transfer_shuffle_heatmap(heatToPlot,
  gg_seuObj = Gg_adult_exc, zf_seuObj = seuObj_glut,
  gg_shuffleThis = "anno_level_3",
  col_use = colorRampPalette(c("white", "#aa7c00", "#6b3e08"))(150),
  anchors = anchors_GGtoZF,
  RDSfiles = RDSfiles,
  rds_name = "lblTrnf_Gluts_lvl3_300shuff.rds"
)

p1 <- res$plot

ggsave(file.path(Figures_dir, "corrplot_glut_lablTrns_annoLvl3_wSig.svg"),
  plot = p1,
  width = 10, height = 9
)

#------ SPEARMAN CORRELATION ------------------------------------------------###
### which genes to use
# genes_fromLblTrans <- anchors_GGtoZF@anchor.features #not used

redoMarkers <- F
if (redoMarkers) {
  Idents(Gg_adult_exc) <- "anno_level_3"
  ggGlut_Markers <- FindAllMarkers(Gg_adult_exc,
    assay = "SCTexc",
    slot = "data",
    logfc.threshold = 0.1,
    min.pct = 0.2,
    min.diff.pct = 0.3,
    only.pos = T,
    recorrect_umi = F
  )
  saveRDS(ggGlut_Markers, paste0(RDSfiles, "/DEGs/ggGlutMarkers_annLvl3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds"))

  Idents(seuObj_glut) <- "clusterName3"
  zfGlut_Markers <- FindAllMarkers(seuObj_glut,
    assay = "SCTexc",
    slot = "data",
    logfc.threshold = 0.1,
    min.pct = 0.2,
    min.diff.pct = 0.3,
    only.pos = T,
    recorrect_umi = F
  )
  saveRDS(zfGlut_Markers, paste0(RDSfiles, "/DEGs/zfGlutMarkers_clusterName3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds"))
}
Idents(Gg_adult_exc) <- "anno_level_3"
Idents(seuObj_glut) <- "clusterName3"

ortho_genes <- intersect(rownames(seuObj_glut), rownames(Gg_adult_exc))

ggGlut_Markers <- readRDS(paste0(RDSfiles, "/DEGs/ggGlutMarkers_annLvl3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds")) %>%
  filter(gene %in% ortho_genes & p_val_adj < 0.05)
zfGlut_Markers <- readRDS(paste0(RDSfiles, "/DEGs/zfGlutMarkers_clusterName3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds")) %>%
  filter(gene %in% ortho_genes & p_val_adj < 0.05)

top_markers_gg <- Extract_Top_Markers(ggGlut_Markers, num_genes = 50, named_vector = FALSE, make_unique = TRUE)
top_markers_zf <- Extract_Top_Markers(zfGlut_Markers, num_genes = 50, named_vector = FALSE, make_unique = TRUE)
top_markers_both <- union(top_markers_gg, top_markers_zf)

# genes_fromWlcDegs <- c(zfGlut_Markers$gene, ggGlut_Markers$gene) %>% unique()
# genes_fromBoth <- intersect(genes_fromLblTrans, genes_fromWlcDegs)

### get the averages
zf_avg <- as.data.frame(AverageExpression(seuObj_glut, assays = "SCTexc", slot = "data", features = top_markers_both))
colnames(zf_avg) <- gsub("SCTexc.", "", colnames(zf_avg))

gg_avg <- as.data.frame(AverageExpression(Gg_adult_exc, assays = "SCTexc", slot = "data", features = top_markers_both))
colnames(gg_avg) <- gsub("SCTexc.", "", colnames(gg_avg))

GSI_heatmap <- generate_SpearmanCorrHeatmap_wSig(
  exprAvg_filt_zf = zf_avg, exprAvg_filt_other = gg_avg,
  zf.levels = glutLvls.byRgnPct,
  other.levels = gg.levels.lvl3.v2,
  minSetTo = 0,
  FigName = paste0("GSI_glut_chkLvl3_varGenesLblTrn"),
  Figures_dir = Figures_dir
)

GSI_heatmap$hm
p2 <- as.ggplot(GSI_heatmap$hm) + coord_equal()
ggsave(file.path(Figures_dir, "corrplot_gsi_annoLvl3.svg"),
  plot = p2,
  width = 10, height = 9
)

#------ COMBINED PLOT -------------------------------------------------------###
colnames(res$sig_label_mat) <- colnames(res$sig_label_mat) %>% gsub("_", ".", .)
rownames(res$sig_label_mat) <- rownames(res$sig_label_mat) %>%
  gsub("-", ".", .) %>%
  gsub("\\(|\\)", ".", .)

comb_mat <- GSI_heatmap$hm@matrix + res$cor_mat
sig_combined_mat <-
  (GSI_heatmap$sig_labels == "*") +
  (res$sig_label_mat == "*")

sig_combined_mat <- ifelse(
  sig_combined_mat == 2, "**",
  ifelse(sig_combined_mat == 1, "*", "")
)


row_max <- apply(comb_mat, 1, max, na.rm = TRUE)

# Redefine cell_fun for combined matrix
star_fun_comb <- function(j, i, x, y, width, height, fill) {
  lab <- sig_combined_mat[i, j]
  if (!is.na(lab) && lab != "") {
    grid.text(lab, x, y, gp = gpar(fontsize = 12, col = "red"))
  }
}

row_anno <- rowAnnotation(
  MaxScore = anno_barplot(row_max,
    gp = gpar(fill = "#555555"),
    width = unit(1.5, "cm")
  )
)

# row_max <- apply(comb_mat, 1, max, na.rm = TRUE)
# comb_mat2 <- cbind(comb_mat, Max = row_max)

# Build the combined heatmap with annotation
comb_heat <- pheatmap(
  comb_mat,
  cell_fun = star_fun_comb,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  col = colorRampPalette(c("white", "#aa7c00", "#6b3e08"))(150),
  right_annotation = row_anno,
  border_color = "grey",
  cellwidth = 15,
  cellheight = 15
)

comb_heat
p3 <- as.ggplot(comb_heat) + coord_equal()

ggsave(file.path(Figures_dir, "corrplot_combined_annoLvl3_w2sig.svg"),
  plot = p3,
  width = 10, height = 9
)
