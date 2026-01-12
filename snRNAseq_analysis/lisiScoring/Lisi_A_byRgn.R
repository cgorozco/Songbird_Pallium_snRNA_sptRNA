# LIBRARIES
library(lisi)
library(ggridges)
library(Seurat)
library(tidyverse)
library(ggpubr)

# DIRECTORIES
root_dir <- '/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/15_LISIscoring'
RDSfiles <- paste0(root_dir, "/RDS_files")
Figures_dir <- paste0(root_dir, "/Figures")

# LOAD DATA AND PREP
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")
seuObj <- readRDS(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")


useLisi_Meta <- FetchData(seuObj, vars = c("clusterName3", "RegionAcronym", "BirdID", "orig.ident"))

if (!file.exists(paste0(RDSfiles, "/lisiScore_byRgn_20251120.rds"))){
  pc30_embeds <- as.data.frame(Embeddings(seuObj, reduction = "integrated.rpca_pca.SCTorigIdent3")[, 1:30])
  
  region_lisi <- compute_lisi(pc30_embeds, useLisi_Meta, c("RegionAcronym"))
  region_lisi <- region_lisi %>% rename("lisiScore" = colnames(region_lisi[1]))
  saveRDS(region_lisi, file = paste0(RDSfiles, "/lisiScore_byRgn_20251120.rds"))
}

region_lisi <- readRDS(file.path(RDSfiles, "lisiScore_byRgn_20251120.rds"))



###===== PLOT AND CALCULATE SIGNIFICANCE ACROSS DIFFERENT CELL CLASSES
lisiScores <- bind_cols(useLisi_Meta, region_lisi)
lisiScores$songRgn <- ifelse(lisiScores$RegionAcronym %in% c("HVC", "LMAN", "RA"), "song", "nonSong") 
lisiScores$songRgn <- factor(lisiScores$songRgn, levels = c("song", "nonSong"))
lisiScores$CellClass <- lisiScores$clusterName3

gluts <- seuObj$clusterName3[grepl("^Glut", seuObj$clusterName3)] %>% unique()
gabas <- seuObj$clusterName3[grepl("MGE|CGE|LGE", seuObj$clusterName3)] %>% unique()
npc <- seuObj$clusterName3[grepl("^Neurogenic", seuObj$clusterName3)] %>% unique()
oligos <-  seuObj$clusterName3[grepl("Oligo|COP|OPC", seuObj$clusterName3)] %>% unique()
epen <-  seuObj$clusterName3[grepl("14|17", seuObj$clusterName3)] %>% unique()

lisiScores <- lisiScores %>% 
  mutate(
    CellClass = case_when(
      CellClass %in% gluts ~ "Exc",
      CellClass %in% gabas ~ "Inh",
      CellClass %in% npc ~ "NPCs",
      CellClass %in% oligos ~ "Olig",
      CellClass %in% epen ~ "Epen",
      TRUE ~ CellClass
    )
  )
lisiScores$CellClass <- factor(lisiScores$CellClass, 
                               levels = c("Exc", "Inh", "NPCs", "Astro(12)", "Epen", "Olig", "Micro(16)", "Endo(15)"))


# visualize
medians_df <- lisiScores %>%
  group_by(CellClass, songRgn) %>%
  summarize(median_lisi = median(lisiScore), .groups = "drop")

ggplot(lisiScores, aes(x = CellClass, y = lisiScore, fill = songRgn)) +
  gghalves::geom_half_violin(
    data = subset(lisiScores, songRgn == "song"),
    aes(x = CellClass),
    side = "l",
    trim = FALSE,
    draw_quantiles = T
  ) +
  gghalves::geom_half_violin(
    data = subset(lisiScores, songRgn == "nonSong"),
    aes(x = CellClass),
    side = "r",
    trim = FALSE
  ) +
  scale_fill_manual(values = c("song" = "#800000", "nonSong" = "gray")) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    axis.line = element_line(color = "#060606"),
    axis.ticks = element_line(color = "#060606")
  ) + 
  NoLegend() +
  scale_y_continuous(breaks = seq(0, ceiling(max(lisiScores$lisiScore)), by = 1)) +
  # Add median value text
  geom_text(
    data = medians_df,
    aes(
      x = CellClass,
      y = median_lisi,
      label = round(median_lisi, 1),
      group = songRgn
    ),
    position = position_nudge(x = ifelse(medians_df$songRgn == "song", -0.3, 0.3)),
    size = 3,
    fontface = "bold"
  ) +
  RotatedAxis()


# stats
lisiScores_avg <- lisiScores %>%
  group_by(CellClass, BirdID, RegionAcronym, songRgn) %>%
  summarise(avg_lisi = mean(lisiScore, na.rm = F), .groups = "drop")

ggboxplot(lisiScores_avg, x = "songRgn", y = "avg_lisi",
          color = "songRgn", palette = "uchicago",
          add = "jitter",
          facet.by = "CellClass", short.panel.labs = FALSE,
          ncol = 9) + stat_compare_means(label = "p.signif")


