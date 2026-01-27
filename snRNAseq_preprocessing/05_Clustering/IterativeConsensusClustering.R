# This uses the scrattch.hicat package developed by the Allen Institute.
# Used to perform iterative clustering several times (50) with 80% of random subsampled cells,
# to produce a consensus higher-confidence clusters.

## ===============================
## Libraries
## ===============================
library(Seurat)
library(tidyverse)
library(Matrix)
library(scrattch.hicat)

## ===============================
## Paths
## ===============================
root_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering"
RDSfiles <- file.path(root_dir, "RDS_files")
dir.create(RDSfiles, showWarnings = FALSE, recursive = TRUE)

## ===============================
## Load Seurat object
## ===============================
seuObj <- readRDS(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/05c_ClassifyCells_subcluster_083024/RDS_files/OutputObject_toUse/snRNAseq_allCells_20241214.rds"
)

## ===============================
## Prepare normalized matrix for HiCat
## ===============================
count_matrix <- GetAssayData(
  seuObj,
  assay = "SCTorigIdent3",
  layer = "counts"
)

norm.dat <- Matrix(cpm(count_matrix), sparse = TRUE)
norm.dat@x <- log2(norm.dat@x + 1)

all.cells <- colnames(norm.dat)

## ===============================
## HiCat parameters
## ===============================
de.param <- de_param(
  padj.th = 0.05,
  lfc.th = 1,
  low.th = 1,
  q1.th = 0.3,
  q.diff.th = 0.7,
  de.score.th = 150
)

## ===============================
## Run consensus clustering
## ===============================
run_consensus_clust(
  norm.dat,
  niter = 50,
  de.param = de.param,
  dim.method = "PCA",
  mc.cores = 4
)

## ===============================
## Collect subsample results
## ===============================
result_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_Spatial/RDS_files/Iter30"
result.files <- file.path(result_dir, dir(result_dir, pattern = "result.*.rda"))

co.result <- collect_subsample_cl_matrix(
  norm.dat,
  result.files,
  all.cells = all.cells
)

## ===============================
## Iterative consensus clustering
## ===============================
consensus.result <- iter_consensus_clust(
  cl.list = co.result$cl.list,
  cl.mat = co.result$cl.mat,
  norm.dat = norm.dat,
  select.cells = all.cells,
  de.param = de.param,
  merge.type = "directional",
  method = "auto"
)

## ===============================
## Save marker genes
## ===============================
markerGenes <- consensus.result[["markers"]]
saveRDS(markerGenes, file = file.path(RDSfiles, "markersIterClust.rds"))

## ===============================
## Refine clusters and add to Seurat
## ===============================
refine.result <- refine_cl(
  consensus.result$cl,
  cl.mat = co.result$cl.mat,
  tol.th = 0.01,
  confusion.th = 0.8,
  min.cells = de.param$min.cells
)

seuObj <- AddMetaData(
  seuObj,
  metadata = as.data.frame(refine.result$cl),
  col.name = "IterConsenClusters"
)

## ===============================
## Build cluster naming table
## ===============================
meta4Lbl <- FetchData(
  seuObj,
  vars = c(
    "IterConsenClusters",
    "cellClass_lvl2",
    "cellClass_lvl3",
    "RegionAcronym",
    "PallRgn"
  )
)

meta4Lbl$PallRgnAlt <- ifelse(
  meta4Lbl$PallRgn %in% c("NidoPall", "HyperPall"),
  "NidoHyperPall",
  as.character(meta4Lbl$PallRgn)
)

pct_thresh <- 80
summary_tables <- list()

for (col in colnames(meta4Lbl)[-1]) {
  pct_tbl <- meta4Lbl %>%
    group_by(IterConsenClusters, !!sym(col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(IterConsenClusters) %>%
    mutate(pct = n / sum(n) * 100)

  summary_tables[[col]] <- pct_tbl %>%
    group_by(IterConsenClusters) %>%
    summarise(
      value = ifelse(
        any(pct >= pct_thresh),
        as.character((!!sym(col))[which.max(pct)]),
        paste0("NotPct", pct_thresh)
      ),
      .groups = "drop"
    ) %>%
    rename(!!paste0("ALN_", col) := value)
}

mergeTable <- reduce(summary_tables, full_join, by = "IterConsenClusters")

mergeTable$ALN_PallRgnUse <- ifelse(
  mergeTable$ALN_PallRgn == paste0("NotPct", pct_thresh),
  mergeTable$ALN_PallRgnAlt,
  mergeTable$ALN_PallRgn
)

mergeTable$clusterName <- paste(
  ifelse(
    mergeTable$ALN_cellClass_lvl3 != paste0("NotPct", pct_thresh),
    mergeTable$ALN_cellClass_lvl3,
    mergeTable$ALN_cellClass_lvl2
  ),
  ifelse(
    mergeTable$ALN_RegionAcronym != paste0("NotPct", pct_thresh),
    mergeTable$ALN_RegionAcronym,
    mergeTable$ALN_PallRgnUse
  ),
  sep = "-"
)

mergeTable$clusterName <- paste0(
  mergeTable$clusterName,
  "(",
  mergeTable$IterConsenClusters,
  ")"
)

## ===============================
## Manual renaming → clusterName3
## ===============================
metaRename <- mergeTable %>%
  select(IterConsenClusters, clusterName) %>%
  mutate(
    clusterName2 = case_when(
      IterConsenClusters == 11 ~ "CGE",
      IterConsenClusters %in% c(47, 48, 49) ~ "CGE_LHX8",
      IterConsenClusters == 3 ~ "LGE",
      IterConsenClusters == 35 ~ "LGE_Pre",
      IterConsenClusters %in% c(50, 51, 52) ~ "MGE_ST18",
      IterConsenClusters %in% c(32, 36, 37, 38) ~ "MGE_SST",
      IterConsenClusters %in% c(33, 39) ~ "MGE_song",
      IterConsenClusters %in% c(40, 31, 34) ~ "MGE_PVALB",
      TRUE ~ clusterName
    ),
    clusterName3 = case_when(
      IterConsenClusters %in% c(58, 55, 54) ~ "COP-NFOL",
      TRUE ~ clusterName2
    )
  )

metaRename$clusterName3 <- paste0(
  metaRename$clusterName3,
  "(",
  metaRename$IterConsenClusters,
  ")"
)

## ===============================
## Add clusterName3 to Seurat
## ===============================
seuObj <- AddMetaData(
  seuObj,
  metadata = metaRename["clusterName3"]
)

## ===============================
## Save final object
## ===============================
saveRDS(
  seuObj,
  file = file.path(RDSfiles, "OutputUse/seuObj_iterCluster_clean_20250131.rds")
)
