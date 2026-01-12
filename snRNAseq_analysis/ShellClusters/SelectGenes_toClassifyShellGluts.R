# LOAD LIBRARIES
library(Seurat)
library(tidyverse)
library(scCustomize)
library(patchwork)

source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")

# LOAD DIRECTORIES
root_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/23_ShellClusters"

RDSfiles <- paste0(root_dir, "/RDS_files")
dir.create(RDSfiles, showWarnings = F)

Figures_dir <- paste0(root_dir, "/Figures")
dir.create(Figures_dir, showWarnings = F)

Tables_dir <- paste0(root_dir, "/Tables")
dir.create(Tables_dir, showWarnings = F)

# LOAD DATA AND PREP
seuObj <- readRDS(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds"
)
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")

matGluts <- coarse_lvl$Cluster[coarse_lvl$Class == "Glutamatergic" & coarse_lvl$Subclass != "Neurogenic"]

seuObj$clusterName3 <- factor(seuObj$clusterName3, levels = glutLvls_byRgnPct)
Idents(seuObj) <- "clusterName3"

group_pcts <- seuObj@meta.data %>%
  dplyr::count(clusterName3, RegionAcronym, name = "n_cells") %>%
  group_by(clusterName3) %>%
  mutate(perc = 100 * n_cells / sum(n_cells)) %>%
  ungroup() # %>% filter(n_cells < 50)
group_pcts$clst_rgn <- paste(group_pcts$clusterName3, group_pcts$RegionAcronym, sep = "_")
passPct <- group_pcts$clst_rgn[group_pcts$perc > 2]


### get only cells that contribute >2% to cluster for clarity
seuObj$clst_rgn <- paste(seuObj$clusterName3, seuObj$RegionAcronym, sep = "_")


## ---- LMAN and AN
seuObj_lman <- subset(seuObj, subset = clusterName3 %in% c(matGluts) & RegionAcronym %in% c("LMAN", "An") & clst_rgn %in% passPct)

table(seuObj_lman$clusterName3)
useThese <- names(which(table(seuObj_lman$clusterName3) > 200))
seuObj_lman_sub <- subset(seuObj_lman, subset = clusterName3 %in% c(useThese))

Idents(seuObj_lman_sub) <- factor(seuObj_lman_sub$clusterName3,
  levels = c(
    "Glut-LMAN(45)", "Glut-LMAN(44)",
    "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(20)",
    "Glut(21)", "Glut-NidoHyperPall(25)"
  )
)

LMANgenesUse <- c(
  "CALCB", "ADSS1", "CNN3",
  "SEMA3E",
  "CAMK2A", "UNC5D"
)

lmanDot <- DotPlot(seuObj_lman_sub,
  cols = c("#800000", "steelblue"), features = LMANgenesUse,
  split.by = "RegionAcronym",
  col.min = 0, col.max = 3, dot.min = 0, dot.max = 100, dot.scale = 10
) + RotatedAxis()

ggsave(file.path(Figures_dir, "dotplot_LMAN_shell.svg"),
  plot = lmanDot,
  height = 7, width = length(LMANgenesUse) * 1.1
)


Stacked_VlnPlot(seuObj_lman_sub,
  features = LMANgenesUse,
  split.by = "RegionAcronym",
  x_lab_rotate = T
)

## ---- RA and Arco
seuObj_ra <- subset(seuObj, subset = clusterName3 %in% c(matGluts) & RegionAcronym %in% c("RA", "Arco") & clst_rgn %in% passPct)

table(seuObj_ra$clusterName3)
useThese <- names(which(table(seuObj_ra$clusterName3) > 300))
seuObj_ra_sub <- subset(seuObj_ra, subset = clusterName3 %in% c(useThese))

seuObj_ra_sub <- PrepSCTFindMarkers(seuObj_ra_sub, assay = "SCTorigIdent3")
arcoMrks <- FindAllMarkers(seuObj_ra_sub,
  min.pct = 0.2, min.diff.pct = 0.3,
  logfc.threshold = 1
)
arcoMrks <- arcoMrks %>% Add_Pct_Diff()

Idents(seuObj_ra_sub) <- factor(seuObj_ra_sub$clusterName3, levels = glutLvls_byRgnPct)
RAgenesUse <- c(
  "SRD5A2", "SLC4A11",
  # "TSHZ1", "TSHZ3",
  # "UNC5A", "ADCYAP1",
  "TENM3", "CNTN5",
  # "CERKL",
  "CACNB2"
)

raDot <- DotPlot(seuObj_ra_sub,
  cols = c("#800000", "steelblue"), features = RAgenesUse,
  split.by = "RegionAcronym",
  col.min = 0, col.max = 3, dot.min = 0, dot.scale = 10
) + RotatedAxis()


ggsave(file.path(Figures_dir, "dotplot_RA_shell.svg"),
  plot = raDot,
  height = 7, width = length(RAgenesUse) * 1.1
)

Stacked_VlnPlot(seuObj_ra_sub,
  features = RAgenesUse,
  split.by = "RegionAcronym",
  x_lab_rotate = T,
  pt.size = 0.1, raster = F
)

## ---- HVC and Nidos
seuObj_hvc <- subset(seuObj, subset = clusterName3 %in% c(matGluts) & RegionAcronym %in% c("HVC", "NC", "NCM") & clst_rgn %in% passPct)

table(seuObj_hvc$clusterName3)
useThese <- names(which(table(seuObj_hvc$clusterName3) > 300))
seuObj_hvc_sub <- subset(seuObj_hvc, subset = clusterName3 %in% c(useThese))

### nido degs
seuObj_hvc_sub <- PrepSCTFindMarkers(seuObj_hvc_sub, assay = "SCTorigIdent3")
nidoMrks <- FindAllMarkers(seuObj_hvc_sub,
  min.pct = 0.2, min.diff.pct = 0.3,
  logfc.threshold = 1
)
nidoMrks <- nidoMrks %>% Add_Pct_Diff()
nidoMrks_list <- split(nidoMrks, f = nidoMrks$cluster)

hvcS_v28 <- setdiff(
  intersect(nidoMrks_list$`Glut-HVC(7)`$gene, nidoMrks_list$`Glut-HVC(8)`$gene),
  nidoMrks_list$`Glut-NidoHyperPall(28)`$gene
)

### visualize
Idents(seuObj_hvc_sub) <- factor(seuObj_hvc_sub$clusterName3, levels = glutLvls_byRgnPct)
HVCgenesUse <- c(
  "SEMA3E", "ADSS1", "WNT5B",
  "ALDH1A2", "CADPS2", "ZEB2",
  "GRIA4", "CACNA1G"
)

hvcDot <- DotPlot(seuObj_hvc_sub,
  cols = c("#800000", "steelblue", "darkblue"), features = HVCgenesUse,
  split.by = "RegionAcronym",
  col.min = 0, col.max = 3, dot.min = 0, dot.scale = 10
) + RotatedAxis()

ggsave(file.path(Figures_dir, "dotplot_HVC_shell.svg"),
  plot = hvcDot,
  height = 7, width = length(HVCgenesUse) * 1.1
)

Stacked_VlnPlot(seuObj_hvc_sub,
  features = HVCgenesUse,
  split.by = "RegionAcronym",
  x_lab_rotate = T,
  pt.size = 0.1, raster = F
)


### --- saveDots together
combined_plot <- lmanDot + hvcDot + raDot +
  plot_layout(widths = c(length(LMANgenesUse), length(HVCgenesUse), length(RAgenesUse))) &
  theme(
    legend.position = "top",
    axis.line = element_line(color = "#030303"),
    axis.ticks = element_line(color = "#030303")
  )

ggsave(file.path(Figures_dir, "dotplot_comb_shell.svg"),
  plot = combined_plot,
  height = 7, width = 15.5
)
