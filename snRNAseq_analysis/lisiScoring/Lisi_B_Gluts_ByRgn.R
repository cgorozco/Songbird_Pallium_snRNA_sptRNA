# LIBRARIES
library(lisi)
library(ggridges)
library(Seurat)
library(tidyverse)
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")

ogId_cols <- setNames(origIdent_colors, origIdent_lvls)
rgn_cols <- setNames(region_colors, lvl_order)

# DIRECTORIES
root_dir <- '/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/15_LISIscoring'
RDSfiles <- paste0(root_dir, "/RDS_files")
Figures_dir <- paste0(root_dir, "/Figures")

# LOAD DATA AND PREP
seuObj <- readRDS(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")


useLisi_Meta <- FetchData(seuObj, vars = c("clusterName3", "RegionAcronym", "BirdID", "orig.ident"))
pc30_embeds <- as.data.frame(Embeddings(seuObj, reduction = "integrated.rpca_pca.SCTorigIdent3")[, 1:30])

region_lisi <- readRDS(file = paste0(RDSfiles, "/lisiScore_byRgn_20251120.rds"))
lisiScores <- bind_cols(useLisi_Meta, region_lisi)

###===== LISI VISUALIZE

### GLUTAMATERGIC
gluts <- seuObj$clusterName3[grepl("^Glut", seuObj$clusterName3)] %>% unique()

lisiScores_use <- lisiScores %>% filter(clusterName3 %in% gluts)
lisiScores_use$clusterName3 <- factor(lisiScores_use$clusterName3 , levels = glutLvls_byRgnPct)

# w/filtering
group_counts <- seuObj@meta.data %>%
  filter(clusterName3 %in% gluts) %>%
  count(clusterName3, RegionAcronym, name = "n_cells") %>%
  group_by(clusterName3) %>%
  mutate(perc = 100 * n_cells / sum(n_cells)) %>%
  ungroup() # %>% filter(n_cells < 50)


group_counts$rgn_clust <- paste0(group_counts$RegionAcronym, group_counts$clusterName3)
useThese <- group_counts$rgn_clust[group_counts$perc > 2] %>% unique()

lisiScores_use$rgn_clust <- paste0(lisiScores_use$RegionAcronym, lisiScores_use$clusterName3)
lisiScores_use_filt <- lisiScores_use %>% filter(rgn_clust %in% useThese)

lisiScores_avg <- lisiScores_use_filt %>%
  group_by(clusterName3, RegionAcronym, BirdID) %>%
  summarise(avg_lisi = mean(lisiScore, na.rm = F), .groups = "drop")

lisiScores_avg$RegionAcronym <- factor(lisiScores_avg$RegionAcronym, levels = lvl_order)
bxplot_byOgId <- ggplot(lisiScores_avg, aes(x = clusterName3, y = avg_lisi)) +
                        geom_boxplot(outlier.shape = NA) +  # avoid plotting outliers twice
          geom_jitter(aes(color = RegionAcronym), width = 0.2, size = 1.5) +
                        scale_color_manual(values = rgn_cols) + # one dot per rgn_id
                        theme_classic() +
                        labs(x = "Cluster", y = "Average LISI score") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                              axis.line = element_line(color = "#020202")) +
                        coord_flip()

ggsave(file.path(Figures_dir, "lisiBxPlt_gluts_Clusters_byOgId.svg"), plot = bxplot_byOgId,
       width = 5, height = 5)

