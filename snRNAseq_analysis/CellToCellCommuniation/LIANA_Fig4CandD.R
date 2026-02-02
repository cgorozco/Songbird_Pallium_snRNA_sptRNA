## this uses the pv_common_up_spec of sst_common_up_spec ouputs of
## Songbird_Pallium_snRNA_sptRNA/snRNAseq_analysis/DifferentialGeneExpression/GABA/21A_DE_gaba_wIn_andToGG_Fig4.R

library(Seurat)
library(scCustomize)
library(tidyverse)
library(liana)

seuObj <- readRDS("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/OutputUse/seuObj_iterCluster_clean_20250131.rds")
source("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/Functions/addMetaTo_seuObj.R")

if (F) {
  liana_test <- liana_wrap(seuObj)
  saveRDS(liana_test, file.path("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/17_CellChat/RDS_files/LIANA/lianaTest.rds"))
}
liana_test <- readRDS(file.path("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/17_CellChat/RDS_files/LIANA/lianaTest.rds"))

# We can aggregate these results into a tibble with consensus ranks
liana_test_2 <- liana_test %>%
  liana_aggregate()

dplyr::glimpse(liana_test_2)
liana_test_2$aggregate_rank %>% median()

liana_test_sig <- filter(liana_test_2, aggregate_rank < 0.05)


liana_test_sig$source %>% unique()
liana_test_sig$target <- factor(liana_test_sig$target, levels = dendro_order_3)
liana_test_sig$source <- factor(liana_test_sig$source, levels = dendro_order_3)
liana_test_sig$receptor.complex %>% factor()

liana_test_sig$sourceTarget <- paste(liana_test_sig$source, liana_test_sig$target,
  sep = "->"
)
liana_test_sig$ligandReceptor <- paste(liana_test_sig$ligand.complex, liana_test_sig$receptor.complex,
  sep = "-"
)

clusters_use_source <- c(33, 39, 40)
clusters_use_target <- c(45, 7, 8, 19, 12, 33, 39, 40) # switched these as clusters_use_source and clusters_use_target

cluster_regex_source <- paste0(
  "\\((?:", paste(clusters_use_source, collapse = "|"), ")\\)"
)

cluster_regex_target <- paste0(
  "\\((?:", paste(clusters_use_target, collapse = "|"), ")\\)"
)


liana_use <- liana_test_sig[
  grepl(cluster_regex_source, liana_test_sig$source) &
    grepl(cluster_regex_target, liana_test_sig$target),
]

liana_trunc_1 <- liana_use %>%
  filter(aggregate_rank <= 0.05) %>%
  filter(ligand.complex %in% c(pv_common_up_spec, sst_common_up_spec)) # changed between ligand.complex and receptor.complex to generate each corresponding plot

## === plots

p1 <- ggplot(
  liana_trunc_1,
  aes(
    x = ligandReceptor,
    y = sourceTarget,
    size = natmi.edge_specificity,
    color = sca.LRscore
  )
) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "viridis") +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Ligand–Receptor",
    y = "Source → Target",
    size = "NATMI edge specificity",
    color = "SCA LR score"
  )

Fig_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/17_CellChat/Figures/LIANA"
ggsave(file.path(Fig_dir, "LIANAdotplot_MGEsong_CCC_source_v2_20251222.svg"),
  plot = p1,
  width = 7, height = 5.5
) # w=5 for target and =7 for source
