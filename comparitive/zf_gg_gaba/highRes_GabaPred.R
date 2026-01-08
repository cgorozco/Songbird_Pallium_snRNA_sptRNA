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

Figures_dir <- paste0(root_dir, "/Figures/gaba")
dir.create(Figures_dir, showWarnings = F)

Tables_dir <- paste0(root_dir, "/Tables")
dir.create(Tables_dir, showWarnings = F)

# LOAD DATA
seuObj <- readRDS("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")
Gg_adult <- readRDS(file = paste0(RDSfiles, "/ChkZaremba25_wOrthos_sct.rds"))


Gg_adult_inh <- subset(Gg_adult, subset = tg_class == "Inhibitory neurons")
Idents(Gg_adult_inh) <- "anno_level_3"

seuObj_gaba <- subset(seuObj, subset = cn3_class %in% "GABAergic")
Idents(seuObj_gaba) <- "clusterName3"

re_do <- F
if (re_do == T) {
  Gg_adult_inh <- SCTransform(Gg_adult_inh,
    assay = "RNA_ortho", new.assay.name = "SCTinh",
    variable.features.n = 6000
  )
  saveRDS(Gg_adult_inh, file = paste0(RDSfiles, "/gg_inh_SCT_20251119.rds"))

  seuObj_gaba <- SCTransform(seuObj_gaba,
    assay = "RNA", new.assay.name = "SCTinh",
    variable.features.n = 6000
  )
  saveRDS(seuObj_gaba, file = paste0(RDSfiles, "/zf_inh_SCT_20251119.rds"))
} else {
  Gg_adult_inh <- readRDS(file = paste0(RDSfiles, "/gg_inh_SCT_20251119.rds"))
  seuObj_gaba <- readRDS(file = paste0(RDSfiles, "/zf_inh_SCT_20251119.rds"))
}


varFeats_AllCells <- intersect(Gg_adult@assays[["SCT_ortho"]]@var.features, seuObj@assays[["SCTorigIdent3"]]@var.features)
varFeats_GabaCells <- intersect(Gg_adult_inh@assays[["SCTinh"]]@var.features, seuObj_gaba@assays[["SCTinh"]]@var.features)
varFeats_both <- union(varFeats_AllCells, varFeats_GabaCells)

redo2 <- F
if (redo2) {
  anchors_GGtoZF <- FindTransferAnchors(
    reference = Gg_adult_inh, query = seuObj_gaba,
    features = varFeats_both,
    reference.assay = "SCTinh",
    query.assay = "SCTinh",
    normalization.method = "SCT",
    recompute.residuals = F,
    reduction = "cca"
  )

  saveRDS(anchors_GGtoZF, file = paste0(RDSfiles, "/anchors_inh_cca_SCTorthoToOrigIdent3_20251021.rds"))
} else {
  anchors_GGtoZF <- readRDS(file = paste0(RDSfiles, "/anchors_inh_cca_SCTorthoToOrigIdent3_20251021.rds"))
}


predictions_join_ann3 <- TransferData(
  anchorset = anchors_GGtoZF, refdata = Gg_adult_inh$anno_level_3,
  weight.reduction = "cca"
)
seuObj_gaba <- AddMetaData(seuObj_gaba, metadata = predictions_join_ann3)

# VISUALIZE
zf_meta <- FetchData(seuObj_gaba, vars = c("clusterName3"))

metaComb <- full_join(
  rownames_to_column(zf_meta, var = "rowname"),
  rownames_to_column(predictions_join_ann3, var = "rowname"),
  by = "rowname"
)

heatToPlot <- metaComb %>%
  select(-prediction.score.max) %>%
  group_by(clusterName3) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames("clusterName3") %>%
  select(where(~ any(. > 0.00)))

colnames(heatToPlot) <- gsub("prediction.score.", "", colnames(heatToPlot))

ggInh_lvl3 <- c(
  "Inh_MGE_like_ANO1", "Inh_MGE_like_LHX8",
  "Inh_CGE_like_CACNA2D3", "Inh_CGE_like_MAN1A1", "Inh_CGE_like_CACNA1B",
  "Inh_MGE_like_SST_ABI3BP", "Inh_MGE_like_SST_SLC24A2",
  "Inh_MGE_like_SST_PVALB", "Inh_MGE_like_ST18_PVALB", "Inh_MGE_like_ST18",
  "Inh_OB_CPA6", "Inh_OB_FGFRL1",
  "Inh_LGE_like_CHODL", "Inh_LGE_like_FOXP1_DGKH", "Inh_LGE_like_FOXP1_SGCG",
  "Inh_LGE_like_FOXP2", "Inh_LGE_like_STK31", "Inh_LGE_like_Pre"
)

heatToPlot <- heatToPlot[order(factor(rownames(heatToPlot), levels = rev(dendro_order_3))), ]
heatToPlot <- heatToPlot[, order(factor(colnames(heatToPlot), levels = rev(ggInh_lvl3)))]

# Gaba_corr_heatmap <- pheatmap(heatToPlot,
#                               cellwidth = 15, cellheight = 15,
#                               cluster_rows = F, cluster_cols = F,
#                               color = colorRampPalette(c("white", "#aa7c00", "#6b3e08"))(100)) # colorRampPalette(c("white", "#aa7c00"))(100)) / viridis(100, option = "viridis") / rev(viridis(100, option = "mako"))
#
# Gaba_corr_heatmap
# lblTrn_plot <- as.ggplot(Gaba_corr_heatmap) + coord_equal()

res <- run_label_transfer_shuffle_heatmap(heatToPlot,
  gg_seuObj = Gg_adult_inh, zf_seuObj = seuObj_gaba, gg_shuffleThis = "anno_level_3",
  anchors = anchors_GGtoZF,
  RDSfiles = RDSfiles,
  rds_name = "lblTrnf_GABAs_lvl3_300shuff.rds"
)

p1 <- res$plot

ggsave(file.path(Figures_dir, "gaba_lblTrn_wSig.svg"),
  plot = p1,
  width = 8, height = 7
)

#------ SPEARMAN CORRELATION ------------------------------------------------###
### which genes to use
genes_fromLblTrans <- anchors_GGtoZF@anchor.features

redoMarkers <- F
if (redoMarkers) {
  Idents(Gg_adult_inh) <- "anno_level_3"
  ggGaba_Markers <- FindAllMarkers(Gg_adult_inh,
    assay = "SCTinh",
    slot = "data",
    logfc.threshold = 0.1,
    min.pct = 0.2,
    min.diff.pct = 0.3,
    only.pos = T,
    recorrect_umi = F
  )
  saveRDS(ggGaba_Markers, paste0(RDSfiles, "/DEGs/ggGabaMarkers_annLvl3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds"))

  Idents(seuObj_gaba) <- "clusterName3"
  zfGaba_Markers <- FindAllMarkers(seuObj_gaba,
    assay = "SCTinh",
    slot = "data",
    logfc.threshold = 0.1,
    min.pct = 0.2,
    min.diff.pct = 0.3,
    only.pos = T,
    recorrect_umi = F
  )
  saveRDS(zfGaba_Markers, paste0(RDSfiles, "/DEGs/zfGabaMarkers_clusterName3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds"))
} else {
  ggGaba_Markers <- readRDS(paste0(RDSfiles, "/DEGs/ggGabaMarkers_annLvl3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds")) %>% filter(gene %in% ortho_genes)
  zfGaba_Markers <- readRDS(paste0(RDSfiles, "/DEGs/zfGabaMarkers_clusterName3_0.2pct-pctDiff0.3-0.1LogFC_20251027.rds")) %>% filter(gene %in% ortho_genes)
}

Idents(Gg_adult_inh) <- "anno_level_3"
Idents(seuObj_gaba) <- "clusterName3"

ortho_genes <- intersect(rownames(seuObj_gaba), rownames(Gg_adult_inh))

ggGaba_Markers_filt <- ggGaba_Markers %>% filter(p_val_adj < 0.05 & gene %in% ortho_genes)
zfGaba_Markers_filt <- zfGaba_Markers %>% filter(p_val_adj < 0.05 & gene %in% ortho_genes)

top_markers_gg <- Extract_Top_Markers(ggGaba_Markers_filt, num_genes = 50, named_vector = FALSE, make_unique = TRUE)
top_markers_zf <- Extract_Top_Markers(zfGaba_Markers_filt, num_genes = 50, named_vector = FALSE, make_unique = TRUE)
top_markers_both <- union(top_markers_gg, top_markers_zf)

# genes_fromWlcDegs <- c(zfGaba_Markers$gene, ggGaba_Markers$gene) %>% unique()
# genes_fromBoth <- intersect(genes_fromLblTrans, genes_fromWlcDegs)

### get the averages
zf_avg <- as.data.frame(AverageExpression(seuObj_gaba, assays = "SCTinh", slot = "data", features = top_markers_both))
colnames(zf_avg) <- gsub("SCTinh.", "", colnames(zf_avg))

gg_avg <- as.data.frame(AverageExpression(Gg_adult_inh, assays = "SCTinh", slot = "data", features = top_markers_both))
colnames(gg_avg) <- gsub("SCTinh.", "", colnames(gg_avg))

ggInh.lvl3 <- gsub("_", ".", ggInh_lvl3)
GSI_heatmap <- generate_SpearmanCorrHeatmap_wSig(
  exprAvg_filt_zf = zf_avg, exprAvg_filt_other = gg_avg,
  zf.levels = dendro.order.3,
  other.levels = rev(ggInh.lvl3),
  minSetTo = 0,
  FigName = paste0("GSI_Gaba_chkLvl3_varGenesLblTrn"),
  Figures_dir = Figures_dir
)

gsi_plot <- as.ggplot(GSI_heatmap$hm) + coord_equal()

ggsave(file.path(Figures_dir, "gaba_gsi.svg"),
  plot = gsi_plot,
  width = 8, height = 7
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

ggsave(file.path(Figures_dir, "cmb_lblTran_gsi_w2sig.svg"),
  plot = p3,
  width = 8, height = 7
)
