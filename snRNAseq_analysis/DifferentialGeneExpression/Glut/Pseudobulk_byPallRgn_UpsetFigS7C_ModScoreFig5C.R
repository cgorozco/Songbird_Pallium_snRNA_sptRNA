# LOAD LIBRARIES
library(UpSetR)
library(ComplexUpset)
library(scCustomize)
library(tidyverse)

source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/geneLists.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/pairwiseDegs_toUpset.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/moduleScore_plots.R")

# LOAD DIRECTORIES
Figures_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Figures/glut"
Tables_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Tables/glut/pallRgn"
dir.create(Tables_dir, recursive = T)
RDS_files <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/RDSfiles/gluts/pallRgn"
dir.create(RDS_files, recursive = T)


# LOAD DATA AND PREP
seuObj <- readRDS(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds"
)
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")

# GET GLUT NEURONS
cluster_use_list <- seuObj$clusterName3[grepl("Glut", seuObj$clusterName3)] %>%
unique() %>%
.[!grepl("Pre", .)]
seuObj_matNeu <- subset(seuObj, subset = clusterName3 %in% cluster_use_list)
seuObj_matNeu$clusterName3 <- factor(seuObj_matNeu$clusterName3, levels = glutLvls_byRgnPct)


# FILTER REPLICATES WITH TOO FEW CELLS, BY REGION
# group_counts <- seuObj_matNeu@meta.data %>%
#   count(clusterName3, orig.ident, name = "n_cells") %>%  filter(n_cells < 50)
# hist(group_counts$n_cells) # to set threshold

group_pcts <- seuObj_matNeu@meta.data %>%
dplyr::count(clusterName3, RegionAcronym, name = "n_cells") %>%
group_by(clusterName3) %>%
mutate(perc = 100 * n_cells / sum(n_cells)) %>%
ungroup() # %>% filter(n_cells < 50)


group_pcts$grouping <- paste0(group_pcts$clusterName3, "-", group_pcts$RegionAcronym)
keep_groups <- group_pcts$grouping[group_pcts$perc > 2]

seuObj_matNeu$grouping <- paste0(seuObj_matNeu$clusterName3, "-", seuObj_matNeu$RegionAcronym)
seuObj_matNeu_filt <- subset(seuObj_matNeu, subset = grouping %in% keep_groups)

seuObj_matNeu_filt$pallialRgn <- case_when(
seuObj_matNeu_filt$RegionAcronym %in% c("NC", "NCM", "An") ~ "Nido",
seuObj_matNeu_filt$RegionAcronym %in% c("Arco") ~ "Arco",
seuObj_matNeu_filt$RegionAcronym %in% c("Meso", "Av") ~ "Meso",
seuObj_matNeu_filt$RegionAcronym == "Hyper" ~ "Hyper",
seuObj_matNeu_filt$RegionAcronym %in% c("HVC", "LMAN", "RA") ~ "Song",
TRUE ~ NA_character_ # fallback for any unmatched values
)
seuObj_matNeu_filt$pallialRgn <- factor(seuObj_matNeu_filt$pallialRgn,
levels = c("Hyper", "Meso", "Nido", "Arco", "Song")
)

## PSEUDOBULK
seuObj_matNeu_filt_use <- subset(seuObj_matNeu_filt, subset = pallialRgn != "Song")
seuObj_matNeu_agg <- AggregateExpression(seuObj_matNeu_filt_use,
group.by = c("pallialRgn", "BirdID"),
return.seurat = TRUE
)
Idents(seuObj_matNeu_agg) <- "pallialRgn"


## DEGS
pall_fields <- seuObj_matNeu_agg$pallialRgn %>% unique()

re_do_mrks_1 <- F
if (re_do_mrks_1) {
  pallial_degs_list <- list()
  for (cluster1 in pall_fields) {
    for (cluster2 in pall_fields) {
      if (cluster1 == cluster2) next

      comp_name <- paste0(cluster1, "_vs_", cluster2)

      message("Running: ", comp_name)

      res <- FindMarkers(
        object   = seuObj_matNeu_agg,
        assay    = "RNA",
        ident.1  = cluster1,
        ident.2  = cluster2,
        test.use = "DESeq2",
        only.pos = F,
        min.pct  = 1
      ) %>%
        tibble::rownames_to_column("gene") %>%
        mutate(
          cluster_1 = cluster1,
          cluster_2 = cluster2
        ) %>%
        filter(p_val_adj < 0.05)

      pallial_degs_list[[comp_name]] <- res
    }
  }
  saveRDS(pallial_degs_list, file = file.path(RDS_files, "DEGs_bulkByPallRgn_excludesSong"))
}
pallial_degs_list <- readRDS(file = file.path(RDS_files, "DEGs_bulkByPallRgn_excludesSong"))

pallial_degs_df <- bind_rows(pallial_degs_list, .id = "comparison") #--save output
# write.csv(pallial_degs_df,
#   file = file.path(Tables_dir, "pallRgn_DEGs_excludesSong.csv"),
#   row.names = FALSE
# )

pallial_degs_list_sig_up <- lapply(pallial_degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & avg_log2FC > 0.2)
})
pallial_pair_up <- prepare_upset_from_pairwise(pallial_degs_list_sig_up, remove_loc = FALSE)
colnames(pallial_pair_up$upset_df) <- factor(colnames(pallial_pair_up$upset_df), levels = c("genes", "Hyper", "Meso", "Nido", "Arco"))
# write.csv(pallial_pair_up$upset_df, file = file.path(Tables_dir, "pallRgn_DEGs_upset_up.csv"), row.names = FALSE)

pallial_degs_list_sig_down <- lapply(pallial_degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & avg_log2FC > 0.2)
})
pallial_pair_down <- prepare_upset_from_pairwise(pallial_degs_list_sig_down, remove_loc = FALSE)
colnames(pallial_pair_down$upset_df) <- factor(colnames(pallial_pair_down$upset_df), levels = c("genes", "Hyper", "Meso", "Nido", "Arco"))
# write.csv(pallial_pair_down$upset_df, file = file.path(Tables_dir, "pallRgn_DEGs_upset_down.csv"), row.names = FALSE)


## VISUALIZE
upsetR_pall <- upset(
  pallial_pair_up$upset_df,
  intersect = pall_fields,
  sort_sets = FALSE,
  sort_intersections = "descending",
  min_size = 1
)

ggsave(file.path(Figures_dir, "upset_byPallField_4PallRgns_wOut_SongInDEG_up.svg"),
  plot = upsetR_pall,
  width = 4, height = 3
)

upsetR_pall_down <- upset(
  pallial_pair_down$upset_df,
  intersect = pall_fields,
  sort_sets = FALSE,
  sort_intersections = "descending",
  min_size = 1
)

ggsave(file.path(Figures_dir, "upset_byPallField_4PallRgns_wOut_SongInDEG_down.svg"),
  plot = upsetR_pall_down,
  width = 4, height = 3
)


### ===== MODULE SCORING BASED ON PSEUDODEGS =================================###
# unique degs
genes_arco_only <- pallial_pair_up$upset_df$genes[
  pallial_pair_up$upset_df$Arco & rowSums(pallial_pair_up$upset_df[, !(names(pallial_pair_up$upset_df) %in% c("genes", "Arco"))]) == 0
]

genes_nido_only <- pallial_pair_up$upset_df$genes[
  pallial_pair_up$upset_df$Nido & rowSums(pallial_pair_up$upset_df[, !(names(pallial_pair_up$upset_df) %in% c("genes", "Nido"))]) == 0
]

genes_meso_only <- pallial_pair_up$upset_df$genes[
  pallial_pair_up$upset_df$Meso & rowSums(pallial_pair_up$upset_df[, !(names(pallial_pair_up$upset_df) %in% c("genes", "Meso"))]) == 0
]

genes_only_list <- list(
  Arco = genes_arco_only,
  Nido = genes_nido_only,
  Meso = genes_meso_only
)

for (region in names(genes_only_list)) {
  seuObj_matNeu_filt <- AddModuleScore(
    object = seuObj_matNeu_filt,
    features = list(genes_only_list[[region]]),
    name = paste0(region, "_only")
  )
}

# # union degs
# for (region in names(pallial_pair_up$genes_union)) {
#   seuObj_matNeu_filt <- AddModuleScore(
#     object = seuObj_matNeu_filt,
#     features = list(pallial_pair_up$genes_union[[region]]),
#     name = region
#   )
# }

# FeaturePlot_scCustom(seuObj_matNeu_filt, features = "Meso1")
seuObj_matNeu_filt$clusterName3 <- as.character(seuObj_matNeu_filt$clusterName3)

mesoClusters <- seuObj_matNeu_filt$clusterName3[grepl("Meso", seuObj_matNeu_filt$clusterName3)] %>% unique()

nidoSongShell <- c("Glut-LMAN(45)", "Glut-LMAN(44)", "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)")
nidoClusters <- seuObj_matNeu_filt$clusterName3[grep("23|20|21|24|26", seuObj_matNeu_filt$clusterName3)] %>% unique()

p1 <- plot_module_scores(seuObj_matNeu_filt,
  mod1 = "Meso_only", mod2 = "Nido_only",
  filterClusters = c(nidoClusters, mesoClusters, nidoSongShell)
)

ggsave(file.path(Figures_dir, "moduleScore_mesoVnido.svg"),
  plot = p1,
  width = 5, height = 5
)


arcoClusters <- c(
  "Glut-Arco(1)", "Glut-RA(18)", "Glut-Arco(10)",
  "Glut-ArcoPall(27)", "Glut-RA(19)"
)

p2 <- plot_module_scores(seuObj_matNeu_filt,
  mod1 = "Meso", mod2 = "Arco",
  filterClusters = c(mesoClusters, arcoClusters)
)
ggsave(file.path(Figures_dir, "moduleScore_mesoVarco.svg"),
  plot = p2,
  width = 5, height = 5
)


### IN THREE DIMENSION (NIDO, ARCO, MESO) AND ALL GLUT CLUSTERS
library(scatterplot3d)
library(plotly)

# union degs
for (region in names(pallial_pair_up$genes_union)) {
  seuObj_matNeu_filt <- AddModuleScore(
    object = seuObj_matNeu_filt,
    features = list(pallial_pair_up$genes_union[[region]]),
    name = region
  )
}

seuObj_matNeu_filt$clusterName3 <- as.character(seuObj_matNeu_filt$clusterName3)

mesoClusters <- seuObj_matNeu_filt$clusterName3[grepl("Meso", seuObj_matNeu_filt$clusterName3)] %>% unique()
nidoSong <- c("Glut-LMAN(45)", "Glut-HVC(7)", "Glut-HVC(8)")
nidoShell <- c("Glut-LMAN(44)", "Glut-NidoHyperPall(28)")
nidoClusters <- seuObj_matNeu_filt$clusterName3[grep("23|20|21|24|26|25", seuObj_matNeu_filt$clusterName3)] %>% unique()

arcoClusters <- c(
  "Glut-Arco(1)", "Glut-RA(18)", "Glut-Arco(10)",
  "Glut-ArcoPall(27)", "Glut-RA(19)"
)

mod1 <- "Meso_only"
mod2 <- "Nido_only"
mod3 <- "Arco_only"

# get module scores, then avg by bird and cluster, and then by cluster
allMods <- FetchData(seuObj_matNeu_filt, vars = c(
  "clusterName3", "BirdID",
  paste0(mod1, "1"),
  paste0(mod2, "1"),
  paste0(mod3, "1")
))
allMods_avg <- allMods %>%
  group_by(clusterName3, BirdID) %>%
  mutate(
    Mod1Avg_byBird = mean(.data[[paste0(mod1, "1")]]),
    Mod2Avg_byBird = mean(.data[[paste0(mod2, "1")]]),
    Mod3Avg_byBird = mean(.data[[paste0(mod3, "1")]])
  ) %>%
  ungroup() %>%
  group_by(clusterName3) %>%
  mutate(
    Mod1Avg_byCluster = mean(Mod1Avg_byBird),
    Mod2Avg_byCluster = mean(Mod2Avg_byBird),
    Mod3Avg_byCluster = mean(Mod3Avg_byBird)
  ) %>%
  ungroup() %>%
  dplyr::select(clusterName3, Mod1Avg_byCluster, Mod2Avg_byCluster, Mod3Avg_byCluster) %>%
  unique()


## plots
allMods_avg$group <- "Other"

# Meso
allMods_avg$group[allMods_avg$clusterName3 %in% mesoClusters] <- "Meso"

# Nido subtypes (more specific first)
allMods_avg$group[allMods_avg$clusterName3 %in% nidoSong] <- "NidoSong"
allMods_avg$group[allMods_avg$clusterName3 %in% nidoShell] <- "NidoShell"

# General Nido (only if not already assigned)
allMods_avg$group[
  allMods_avg$clusterName3 %in% nidoClusters &
    allMods_avg$group == "Other"
] <- "Nido"

# Arco
allMods_avg$group[allMods_avg$clusterName3 %in% arcoClusters] <- "Arco"

# cols <- c(
#   Meso      = "#7A5E47",
#   NidoSong  = "#830909",
#   NidoShell = "#E8A6B0",
#   Nido      = "#90ADBA",
#   Arco      = "#C35139",
#   Other     = "gray80"
# )

cols <- c(
  Meso      = "#58D658",
  NidoSong  = "#D7191C",
  NidoShell = "#F7C6D0",
  Nido      = "#2C7FB8",
  Arco      = "#FF7F00",
  Other     = "#7A7A7A"
)


p <- plot_ly(
  allMods_avg,
  x = ~Mod1Avg_byCluster,
  y = ~Mod2Avg_byCluster,
  z = ~Mod3Avg_byCluster,
  color = ~group,
  colors = cols,
  text = ~clusterName3,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 9) # line = list(width = 0.5, color = "black")
) %>%
  layout(
    title = list(
      text = "Pallial Module Space",
      font = list(size = 22)
    ),
    scene = list(
      xaxis = list(
        title = list(text = paste0(mod1, " module"), font = list(size = 18)),
        tickfont = list(size = 14)
      ),
      yaxis = list(
        title = list(text = paste0(mod2, " module"), font = list(size = 18)),
        tickfont = list(size = 14)
      ),
      zaxis = list(
        title = list(text = paste0(mod3, " module"), font = list(size = 18)),
        tickfont = list(size = 14)
      )
    )
  )


htmlwidgets::saveWidget(p, file.path(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Figures/glut",
  "rotating_module_plot_uniqueDEGs_v1.html"
))
