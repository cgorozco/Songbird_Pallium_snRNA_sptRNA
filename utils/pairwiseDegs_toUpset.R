prepare_upset_from_pairwise <- function(
    degs_list,
    remove_loc = TRUE,
    cluster_pattern_exclude = NULL
) {
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  
  ##---------------------------------------------------------
  ## 1. Get cluster1 from names (before "_vs_")
  ##---------------------------------------------------------
  cluster1_names <- str_extract(names(degs_list), "^[^_]+")
  degs_by_cluster1 <- split(degs_list, cluster1_names)
  
  ##---------------------------------------------------------
  ## 2. Extract list of gene vectors for each cluster1
  ##---------------------------------------------------------
  genes_by_cluster1 <- lapply(degs_by_cluster1, function(df_list) {
    lapply(df_list, `[[`, "gene")
  })
  
  ##---------------------------------------------------------
  ## 3. Union and intersection for each cluster1
  ##---------------------------------------------------------
  genes_intersect <- lapply(genes_by_cluster1, function(gene_lists) {
    Reduce(intersect, gene_lists)
  })
  
  genes_union <- lapply(genes_by_cluster1, function(gene_lists) {
    Reduce(union, gene_lists)
  })
  
  ##---------------------------------------------------------
  ## 4. Optional filtering:
  ##    A) remove unwanted clusters from names
  ##    B) remove LOC genes
  ##---------------------------------------------------------
  genes_union_use <- genes_union
  
  if (!is.null(cluster_pattern_exclude)) {
    genes_union_use <- genes_union_use[!grepl(cluster_pattern_exclude, names(genes_union_use))]
  }
  
  if (remove_loc) {
    genes_union_use <- lapply(genes_union_use, function(glist) {
      glist[!grepl("^LOC", glist)]
    })
  }
  
  ##---------------------------------------------------------
  ## 5. Convert to UpSet-ready binary matrix
  ##---------------------------------------------------------
  genes_union_df <- genes_union_use %>%
    enframe(name = "cluster", value = "genes") %>%
    unnest(genes) %>%
    distinct() %>%
    mutate(present = TRUE) %>%
    pivot_wider(
      names_from = cluster,
      values_from = present,
      values_fill = FALSE
    )
  
  ## Return everything useful
  list(
    genes_union = genes_union_use,
    genes_intersect = genes_intersect,
    upset_df = genes_union_df
  )
}
