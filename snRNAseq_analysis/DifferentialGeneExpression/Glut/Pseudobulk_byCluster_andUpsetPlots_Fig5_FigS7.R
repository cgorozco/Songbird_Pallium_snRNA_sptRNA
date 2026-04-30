# LOAD LIBRARIES
library(Seurat)
library(UpSetR)
library(ComplexUpset)
library(scCustomize)
library(tidyverse)

source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/CommonLevelsColors.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/geneLists.R")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/pairwiseDegs_toUpset.R")

# LOAD DIRECTORIES
Figures_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Figures/glut"
Tables_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Tables/glut"
dir.create(Tables_dir, recursive = T)
RDS_files <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/RDSfiles/gluts"
dir.create(RDS_files, recursive = T)
date_stamp <- format(Sys.Date(), "%Y-%m-%d")

# LOAD DATA AND PREP
seuObj <- readRDS(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds"
)
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")

re_do_mrks <- F

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


## PSEUDOBULK
seuObj_matNeu_agg <- AggregateExpression(seuObj_matNeu_filt,
  group.by = c("clusterName3", "orig.ident", "RegionAcronym", "BirdID"),
  return.seurat = TRUE
)
seuObj_matNeu_agg$grouping <- paste0(seuObj_matNeu_agg$clusterName3, "-", seuObj_matNeu_agg$RegionAcronym)

seuObj_matNeu_agg$clusterName3 <- factor(seuObj_matNeu_agg$clusterName3, levels = glutLvls_byRgnPct)
Idents(seuObj_matNeu_agg) <- "clusterName3"


## DEGs; for the upset plots
#---- arco
arco_clusters <- c("Glut-RA(19)", "Glut-RA(18)", "Glut-ArcoPall(27)", "Glut-Arco(1)", "Glut-Arco(10)")

re_do_mrks_1 <- F
if (re_do_mrks_1) {
  arco_degs_list <- list()
  for (cluster1 in arco_clusters) {
    for (cluster2 in arco_clusters) {
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

      arco_degs_list[[comp_name]] <- res
    }
  }
  saveRDS(arco_degs_list, file.path(RDS_files, "arco_AllDegs.rds"))
}

arco_degs_list <- readRDS(file.path(RDS_files, "arco_AllDegs.rds"))

arco_degs_df <- bind_rows(arco_degs_list, .id = "comparison") #--save output
write.csv(arco_degs_df,
  file = file.path(Tables_dir, "arco_DEGs_all_comparisons.csv"),
  row.names = FALSE
)

arco_degs_list_sig_up <- lapply(arco_degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & avg_log2FC > 0)
})

arco_degs_list_sig_down <- lapply(arco_degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & avg_log2FC < 0)
})

both_upNdown <- list(up = arco_degs_list_sig_up, down = arco_degs_list_sig_down)

upset_results_arco <- list()
for (dir_name in names(both_upNdown)) {
  degsUse <- both_upNdown[[dir_name]]
  arco_pair <- prepare_upset_from_pairwise(degsUse, remove_loc = FALSE)
  arco_pair$upset_df <- arco_pair$upset_df %>%
    mutate(
      Ortho = case_when(
        genes %in% chkNzf_orthos ~ "Ortho",
        TRUE ~ "Other"
      )
    )
  
  upset_results_arco[[dir_name]] <- list(arco_pair = arco_pair)

  set_cols <- arco_clusters

  upset_plot <- ComplexUpset::upset(
    arco_pair$upset_df,
    intersect = set_cols,
    min_size = 1,
    sort_sets = FALSE,
    sort_intersections = FALSE,
    intersections = "all"
  )
  upset_results_arco[[dir_name]]$complex_upset <- upset_plot

  upsetR_plot <- upset(
    arco_pair$upset_df,
    intersect = set_cols,
    sort_sets = FALSE,
    sort_intersections = FALSE,
    min_size = 3,
    # set_sizes = upset_set_size(
    #   mapping = aes(fill = Ortho)
    # ) +
    #   geom_bar(position = "fill") +
    #   scale_fill_manual(values = c(
    #     "Ortho" = "orange",
    #     "Other" = "grey70"
    #   )) +
    #  scale_y_continuous(labels = scales::percent),
    annotations = list(
      'Set Size Ortho' = list(
        aes = aes(x = intersection, fill = Ortho),
        geom = list(
          geom_bar(position = "fill"),
          scale_y_continuous(labels = scales::percent),
          scale_fill_manual(values = c(
            "Ortho" = "orange",
            "Other" = "grey70"
          )),
          labs(y = "Percentage")
        )
      )
    ),
    intersections = list(
      "Glut-RA(19)", "Glut-RA(18)", "Glut-ArcoPall(27)", "Glut-Arco(1)", "Glut-Arco(10)",
      c("Glut-RA(19)", "Glut-RA(18)"),
      c("Glut-RA(19)", "Glut-ArcoPall(27)"),
      c("Glut-RA(19)", "Glut-RA(18)", "Glut-ArcoPall(27)"),
      c("Glut-RA(19)", "Glut-RA(18)", "Glut-ArcoPall(27)", "Glut-Arco(1)"),
      c("Glut-RA(18)", "Glut-ArcoPall(27)", "Glut-Arco(1)", "Glut-Arco(10)"),
      c("Glut-ArcoPall(27)", "Glut-Arco(1)", "Glut-Arco(10)"),
      c("Glut-Arco(1)", "Glut-Arco(10)")
    )
  )
  upset_results_arco[[dir_name]]$upsetR <- upsetR_plot
  rm(arco_pair)

  # save df for upset
  write.csv(upset_results_arco[[dir_name]][["arco_pair"]][["upset_df"]],
    file = paste0(Tables_dir, "/arco_upsetDEGs_", dir_name, date_stamp, ".csv"), row.names = FALSE
  )
  saveRDS(upset_results_arco, file.path(RDS_files, paste0("upsetDF_arco", date_stamp, ".rds")))
  }

upset_results_arco <- readRDS(file.path(RDS_files, "upsetDF_arco.rds"))

arco_up <- upset_results_arco$up$upsetR
ggsave(file.path(Figures_dir, "upset_byPallField_arcoUp_wOrtho.svg"),
  plot = arco_up,
  width = 6, height = 3
)

arco_down <- upset_results_arco$down$upsetR
ggsave(file.path(Figures_dir, "upset_byPallField_arcoDown_wOrtho.svg"),
  plot = arco_down,
  width = 6, height = 3
)


#---- nido
nido_clusters <- levels(Idents(seuObj_matNeu_agg)) %>%
  .[!grepl("Arco|RA|Meso|\\(9\\)", .)]


re_do_mrks_2 <- F
if (re_do_mrks_2) {
  nido_degs_list <- list()
  for (cluster1 in nido_clusters) {
    for (cluster2 in nido_clusters) {
      if (cluster1 == cluster2) next

      comp_name <- paste0(cluster1, "_vs_", cluster2)

      message("Running: ", comp_name)

      res <- FindMarkers(
        object   = seuObj_matNeu_agg,
        ident.1  = cluster1,
        assay    = "RNA",
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

      nido_degs_list[[comp_name]] <- res
    }
  }
  saveRDS(nido_degs_list, file.path(RDS_files, "nido_AllDegs.rds"))
}

nido_degs_list <- readRDS(file.path(RDS_files, "nido_AllDegs.rds"))

nido_degs_df <- bind_rows(nido_degs_list, .id = "comparison") #--save output
write.csv(nido_degs_df,
  file = file.path(Tables_dir, "nido_DEGs_all_comparisons.csv"),
  row.names = FALSE
)

nido_degs_list_sig_up <- lapply(nido_degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & avg_log2FC > 0)
})
nido_degs_list_sig_down <- lapply(nido_degs_list, function(df) {
  filter(df, p_val_adj < 0.05 & avg_log2FC < 0)
})

both_upNdown <- list(up = nido_degs_list_sig_up, down = nido_degs_list_sig_down)
upset_results_nido <- list()

for (dir_name in names(both_upNdown)) {
  degsUse <- both_upNdown[[dir_name]]
  nido_pair <- prepare_upset_from_pairwise(degsUse, remove_loc = T)
  nido_pair$upset_df <- nido_pair$upset_df %>%
    mutate(
      Ortho = case_when(
        genes %in% chkNzf_orthos ~ "Ortho",
        TRUE ~ "Other"
      )
    )
  
  upset_results_nido[[dir_name]] <- list(nido_pair = nido_pair)

  set_cols <- nido_clusters

  upsetR_plot <- upset(
    nido_pair$upset_df,
    intersect = set_cols,
    sort_sets = FALSE,
    sort_intersections = FALSE,
    min_size = 1,
    annotations = list(
      'Set Size Ortho' = list(
        aes = aes(x = intersection, fill = Ortho),
        geom = list(
          geom_bar(position = "fill"),
          scale_y_continuous(labels = scales::percent),
          scale_fill_manual(values = c(
            "Ortho" = "orange",
            "Other" = "grey70"
          )),
          labs(y = "Percentage")
        )
      )
    ),
    intersections = list(
      "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)", "Glut-NidoHyperPall(23)",
      "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(26)", "Glut-NidoHyperPall(20)",
      "Glut-LMAN(45)", "Glut-LMAN(44)",
      "Glut(21)", "Glut-NidoHyperPall(25)",
      c("Glut-HVC(7)", "Glut-HVC(8)"),
      c("Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)"),
      # c("Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)", "Glut-NidoHyperPall(23)"),
      c("Glut-LMAN(44)", "Glut-LMAN(45)"),
      c("Glut-HVC(8)", "Glut-LMAN(45)"),
      c("Glut-HVC(7)", "Glut-LMAN(45)"),
      c("Glut-HVC(7)", "Glut-HVC(8)", "Glut-LMAN(45)"),
      c("Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)", "Glut-LMAN(44)", "Glut-LMAN(45)"),
      c("Glut-LMAN(44)", "Glut-NidoHyperPall(28)"),
      # c("Glut-NidoHyperPall(23)", "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(26)", "Glut(21)"),
      # c("Glut-NidoHyperPall(23)", "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(26)", "Glut(21)", "Glut-NidoHyperPall(25)"),
      c(
        "Glut-NidoHyperPall(28)", "Glut-NidoHyperPall(23)",
        "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(26)", "Glut-NidoHyperPall(20)",
        "Glut-LMAN(44)",
        "Glut(21)", "Glut-NidoHyperPall(25)"
      ),
      c(
        "Glut-NidoHyperPall(23)",
        "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(26)", "Glut-NidoHyperPall(20)",
        "Glut(21)", "Glut-NidoHyperPall(25)"
      )
    )
  )
  upset_results_nido[[dir_name]]$upsetR <- upsetR_plot
  rm(nido_pair)

  # save df for upset
  write.csv(upset_results_nido[[dir_name]][["nido_pair"]][["upset_df"]],
    file = paste0(Tables_dir, "/nido_upsetDEGs_", dir_name, ".csv"), row.names = FALSE
  )
  saveRDS(upset_results_nido, file.path(RDS_files, "upsetDF_nido.rds"))
  
}

upset_results_nido <- readRDS(file.path(RDS_files, "upsetDF_nido.rds"))

nido_up <- upset_results_nido$up$upsetR
ggsave(file.path(Figures_dir, "upset_byPallField_nidoUp_wOrtho.svg"),
  plot = nido_up,
  width = 6, height = 5
)

nido_down <- upset_results_nido$down$upsetR
ggsave(file.path(Figures_dir, "upset_byPallField_nidoDown_wOrtho.svg"),
  plot = nido_down,
  width = 6, height = 5
)


#### === another upset with just song/shell clusters to check that overlap
set_cols <- c(
  "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)",
  "Glut-LMAN(44)", "Glut-LMAN(45)",
  "Glut-ArcoPall(27)", "Glut-RA(18)", "Glut-RA(19)"
)

combUse <- list()
for (dir in c("up", "down")) {
  ## ---- select data ----
  nidoUse <- upset_results_nido[[dir]]$nido_pair$upset_df %>%
    select(
      genes,
      "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)",
      "Glut-LMAN(44)", "Glut-LMAN(45)"
    )

  arcoUse <- upset_results_arco[[dir]]$arco_pair$upset_df %>%
    select(
      genes,
      "Glut-ArcoPall(27)", "Glut-RA(18)", "Glut-RA(19)"
    )

  ## ---- combine ----
  combUse[[dir]] <- full_join(nidoUse, arcoUse, by = "genes") %>%
    mutate(across(-genes, ~ replace_na(.x, FALSE))) %>%
    filter(if_any(-genes, identity)) %>%
    mutate(Ortho = case_when(
      genes %in% chkNzf_orthos ~ "Ortho",
      TRUE ~ "Other"
    ))

  write.csv(combUse[[dir]],
    file = paste0(Tables_dir, "/songGluts_sharedDEGs_", dir, ".csv"), row.names = FALSE
  )

  ## ---- upset ----
  upsetShareXPallial <- ComplexUpset::upset(
    combUse[[dir]],
    intersect = set_cols,
    sort_sets = FALSE,
    sort_intersections = FALSE,
  min_size = 1,
    annotations = list(
      'Set Size Ortho' = list(
        aes = aes(x = intersection, fill = Ortho),
        geom = list(
          geom_bar(position = "fill"),
          scale_y_continuous(labels = scales::percent),
          scale_fill_manual(values = c(
            "Ortho" = "orange",
            "Other" = "grey70"
          )),
          labs(y = "Percentage")
        )
      )
    ),
    intersections = list(
      "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)",
      "Glut-LMAN(44)", "Glut-LMAN(45)",
      "Glut-ArcoPall(27)", "Glut-RA(18)", "Glut-RA(19)",

      # shared across pallial regions song
      c("Glut-HVC(7)", "Glut-RA(19)"),
      c("Glut-HVC(8)", "Glut-RA(19)"),
      c("Glut-LMAN(45)", "Glut-RA(19)"),
      c("Glut-HVC(7)", "Glut-HVC(8)", "Glut-RA(19)"),
      c("Glut-HVC(7)", "Glut-LMAN(45)", "Glut-RA(19)"),
      c("Glut-HVC(8)", "Glut-LMAN(45)", "Glut-RA(19)"),
      c("Glut-HVC(7)", "Glut-HVC(8)", "Glut-LMAN(45)", "Glut-RA(19)"),

      # shared across pallial regions song-surround
      c("Glut-NidoHyperPall(28)", "Glut-RA(18)"),
      c("Glut-LMAN(44)", "Glut-RA(18)"),
      c("Glut-NidoHyperPall(28)", "Glut-ArcoPall(27)"),
      c("Glut-LMAN(44)", "Glut-ArcoPall(27)"),
      c("Glut-NidoHyperPall(28)", "Glut-LMAN(44)", "Glut-RA(18)"),
      c("Glut-NidoHyperPall(28)", "Glut-LMAN(44)", "Glut-ArcoPall(27)"),
      c("Glut-NidoHyperPall(28)", "Glut-RA(18)", "Glut-ArcoPall(27)"),
      c("Glut-LMAN(44)", "Glut-RA(18)", "Glut-ArcoPall(27)"),
      c("Glut-NidoHyperPall(28)", "Glut-LMAN(44)", "Glut-RA(18)", "Glut-ArcoPall(27)")
    )
  ) 
  
  # ---- save ----
  ggsave(
    file.path(
      Figures_dir,
      paste0("upset_byPallField_Shared_SongOrSorround_wOrthos_", dir, ".svg")
    ),
    plot = upsetShareXPallial,
    width = 6.5,
    height = 4.5
  )
}
upsetShareXPallial
