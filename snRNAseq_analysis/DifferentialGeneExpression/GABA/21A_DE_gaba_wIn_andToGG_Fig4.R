# ========== SET UP
library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(readr)
library(readxl)
library(EnhancedVolcano)
library(ggrastr)
library(writexl)

source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/moduleScore_plots.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")

options(future.globals.maxSize = 100000 * 1024^10)

# SET UP DIRECTORIES
Fig_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Figures"
Tables_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Tables/gaba"
dir.create(Tables_dir, recursive = T)

# LOAD DATA
seuObj <- readRDS("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")

DimPlot(seuObj, group.by = "Phase", label = T) + NoLegend()


# ========== SUBSET TO NEEDED
### ZF
gabas <- coarse_lvl$Cluster[coarse_lvl$Class == "GABAergic"]
seuObj_gaba <- subset(seuObj, subset = clusterName3 %in% c(gabas))
table(seuObj_gaba$clusterName3)

# filter replicates with too few cells
group_counts <- seuObj_gaba@meta.data %>%
  count(BirdID, clusterName3, cn3_subclass, name = "n_cells") # %>%  filter(n_cells < 50)
hist(group_counts$n_cells) # to set threshold

group_counts$grouping <- paste0(group_counts$BirdID, "-", group_counts$clusterName3, "-", group_counts$cn3_subclass)
keep_groups <- group_counts$grouping[group_counts$n_cells > 20] # histogram is 10, but example in scBestPractices is 30, so going w/20

seuObj_gaba$grouping <- paste0(seuObj_gaba$BirdID, "-", seuObj_gaba$clusterName3, "-", seuObj_gaba$cn3_subclass)
seuObj_gaba_agg_1 <- subset(seuObj_gaba, subset = grouping %in% keep_groups)

seuObj_gaba_agg <- AggregateExpression(seuObj_gaba_agg_1,
  group.by = c("BirdID", "clusterName3", "cn3_subclass"),
  assays = "SCTorigIdent3",
  slot = "counts",
  return.seurat = TRUE
)
Idents(seuObj_gaba_agg) <- "cn3_subclass"

### GG
Gg_adult <- readRDS(file = paste0(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/13_chicken",
  "/RDS_files/ChkZaremba25_wOrthos_sct.rds"
))

gaba_ofI <- Gg_adult@meta.data %>%
  distinct(anno_level_3) %>%
  filter(grepl("SST|PVALB", anno_level_3)) %>%
  pull(anno_level_3)

gg_gabaOfI <- subset(Gg_adult, subset = anno_level_3 %in% c(gaba_ofI))

# drop replicates if too few cells
group_counts_gg <- gg_gabaOfI@meta.data %>%
  count(individual, anno_level_3, name = "n_cells") # %>% filter(n_cells < 80)
# won't filter since would remove all of one cell type

gg_gabaOfI_agg <- AggregateExpression(gg_gabaOfI,
  group.by = c("individual", "anno_level_3"), # "anno_level_3"
  assays = "SCT_ortho",
  slot = "counts",
  return.seurat = T
)

# ========== DGE
## Within ZF
clustsSong <- c("PVALB-song", "MGE-song")
clustersToCompare <- c("PVALB", "SST")
degs_list <- list()
for (cluster in clustsSong) {
  for (cluster2 in clustersToCompare) {
    degs_list[[paste0(cluster, "_vs_", cluster2)]] <- FindMarkers(
      object = seuObj_gaba_agg, assay = "SCTorigIdent3",
      ident.1 = cluster,
      ident.2 = cluster2,
      test.use = "DESeq2"
    ) %>% mutate(
      cluster_1 = cluster,
      cluster_2 = cluster2
    )

    degs_list[[paste0(cluster, "_vs_", cluster2)]]$gene <- rownames(degs_list[[paste0(cluster, "_vs_", cluster2)]])
  }
}

degs_list_sig <- lapply(degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & !grepl("^LOC", gene) & avg_log2FC > 0)
})

# common_genes <- Reduce(intersect, lapply(degs_list_sig, `[[`, "gene"))

## Between zf and gg
gg_gabaOfI_agg <- RenameAssays(gg_gabaOfI_agg, SCT_ortho = "SCT")
seuObj_gaba_agg <- RenameAssays(seuObj_gaba_agg, SCTorigIdent3 = "SCT")

agg_combined <- merge(seuObj_gaba_agg, gg_gabaOfI_agg)
agg_combined <- JoinLayers(agg_combined)

agg_combined$clusterComb <- ifelse(is.na(agg_combined$cn3_subclass),
  "gg_PV_SST", agg_combined$cn3_subclass
)
Idents(agg_combined) <- "clusterComb"

degsToGG_list <- list()
cluster2 <- "gg_PV_SST"
for (cluster in clustsSong) {
  degsToGG_list[[paste0(cluster, "_vs_", cluster2)]] <- FindMarkers(
    object = agg_combined, assay = "SCT",
    # features = "GAD2",
    slot = "counts",
    ident.1 = cluster,
    ident.2 = cluster2,
    test.use = "DESeq2"
  )
  degsToGG_list[[paste0(cluster, "_vs_", cluster2)]]$gene <- rownames(degsToGG_list[[paste0(cluster, "_vs_", cluster2)]])
}

degsToGG_list_sig <- lapply(degsToGG_list, function(df) {
  filter(df, p_val_adj < 0.05 & !grepl("^LOC", gene) & avg_log2FC > 0)
})

degsToGG_ctrl_list <- list()
cluster2 <- "gg_PV_SST"
for (cluster in clustersToCompare) {
  degsToGG_ctrl_list[[paste0(cluster, "_vs_", cluster2)]] <- FindMarkers(
    object = agg_combined, assay = "SCT",
    # features = "GAD2",
    slot = "counts",
    ident.1 = cluster,
    ident.2 = cluster2,
    test.use = "DESeq2"
  )
  degsToGG_ctrl_list[[paste0(cluster, "_vs_", cluster2)]]$gene <- rownames(degsToGG_ctrl_list[[paste0(cluster, "_vs_", cluster2)]])
}

degsToGG_ctrl_list_sig <- lapply(degsToGG_ctrl_list, function(df) {
  filter(df, p_val_adj < 0.05 & !grepl("^LOC", gene)) # this not picking up/down direction,
}) # since I want expression equal between two


### COMBINE LISTS
comb_list <- c(degsToGG_list_sig, degs_list_sig)
diff_nonSong <- Reduce(intersect, lapply(degsToGG_ctrl_list_sig, `[[`, "gene"))

pv_song_list <- comb_list[grepl("PVALB-song", names(comb_list))]
pv_common_up <- Reduce(intersect, lapply(pv_song_list, `[[`, "gene"))
pv_common_up_spec <- setdiff(pv_common_up, diff_nonSong)

sst_song_list <- comb_list[grepl("MGE-song", names(comb_list))]
sst_common_up <- Reduce(intersect, lapply(sst_song_list, `[[`, "gene"))
sst_common_up_spec <- setdiff(sst_common_up, diff_nonSong)

up_both <- intersect(pv_common_up_spec, sst_common_up_spec)

#### check
DotPlot_scCustom(seuObj_gaba,
  features = up_both,
  group.by = "cn3_subclass"
) + RotatedAxis()

# VISUALIZE
other_genes <- c("RUNX1", "CUX1", "SETBP1", "NACC2", "SLIT3", "PLNXA4", "ITPR1", "HS3ST5", "OSBPL6")
plot_list <- list()
excel_list <- list()

for (i in seq_along(degs_list)) {
  df <- degs_list[[i]]
  excel_list[[ names(degs_list)[i] ]] <- df

  # make your custom colors
  point_colors <- point_colors_simple(df)
  genes_to_label <- top_genes_label(df, 3, 3, 2, 2)

  plot_list[[names(degs_list)[i]]] <- EnhancedVolcano(
    df,
    lab = rownames(df),
    selectLab = c(genes_to_label, other_genes), # genes_to_label,
    colCustom = point_colors,
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 0,
    drawConnectors = T,
    arrowheads = FALSE,
    legendPosition = "none",
    gridlines.minor = FALSE,
    gridlines.major = FALSE,
    max.overlaps = 15,
    raster = F, # False for tiff file
    # ylim = c(0, max(-log10(df[["p_val_adj"]]), na.rm = TRUE) + 50),

    pointSize = 0.1,
    cutoffLineCol = "#020201",
    cutoffLineWidth = 0.2,
    title = names(degs_list[i]),
    axisLabSize = 6,
    titleLabSize = 4,
    subtitleLabSize = 0,
    captionLabSize = 0,
    labSize = 2,
    widthConnectors = 0.1
  ) + theme(
    axis.line = element_line(size = 0.4, color = "#060606"),
    axis.ticks = element_line(size = 0.4, color = "#060606")
  ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) # removes buffer at y=0, but keeps at the top

  plot_list[[names(degs_list)[i]]] <- rasterize(plot_list[[names(degs_list)[i]]],
    layers = "Point",
    dpi = 600, dev = "ragg"
  )

  ggsave(paste0(Fig_dir, "/vlc_", names(degs_list)[i], ".svg"),
    plot = plot_list[[names(degs_list)[i]]],
    width = 2.2, height = 3
  )
}

write_xlsx(excel_list, path = paste0(Tables_dir, "/wInZF_songPVnSST_otherPVnSST.xlsx"))

comb_plot <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(Fig_dir, "vlcPlotComb_MGEsong_vOtherMGE.svg"),
  plot = comb_plot,
  width = 4, height = 5.7, dpi = 600
)
