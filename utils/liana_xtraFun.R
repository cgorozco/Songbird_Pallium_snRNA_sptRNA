library(dplyr)
library(purrr)

# Main wrapper function
analyze_liana <- function(liana_df, 
                          focus_cluster,    # e.g. "Astro(12)"
                          focus_as = c("target", "source"),  # flip perspective
                          keep_pattern = "Glut",             # which clusters to keep
                          exclude_clusters = NULL,           # e.g. c("23","28","20")
                          sources_of_interest = NULL         # which sources to compare
) {
  
  focus_as <- match.arg(focus_as)
  
  # Ensure lig_rec_pair exists
  liana_df <- liana_df %>%
    mutate(lig_rec_pair = paste(ligand.complex, receptor.complex, sep = "_"))
  
  # Filter based on source or target depending on `focus_as`
  if (focus_as == "target") {
    liana_use <- liana_df %>% filter(target == focus_cluster)
  } else {
    liana_use <- liana_df %>% filter(source == focus_cluster)
  }
  
  # Build list of lig_rec_pair grouped by the opposite axis (source if focus_as=target, target if focus_as=source)
  group_col <- if (focus_as == "target") "source" else "target"
  
  list_targets_use <- liana_use %>%
    filter(grepl(keep_pattern, !!sym(group_col))) %>%
    group_split(!!sym(group_col)) %>%
    set_names(map_chr(., ~ unique(.x[[group_col]]))) %>%
    map(~ .x$lig_rec_pair)
  
  # Union of exclude clusters (if provided)
  targets_union <- NULL
  if (!is.null(exclude_clusters)) {
    clusters_exclude <- unique(liana_use[[group_col]][grepl(paste(exclude_clusters, collapse = "|"), 
                                                            liana_use[[group_col]])])
    targets_union <- unique(unlist(list_targets_use[clusters_exclude]))
  }
  
  # Compare sources of interest
  if (!is.null(sources_of_interest)) {
    result <- map(sources_of_interest, function(src) {
      vals <- list_targets_use[[src]]
      if (!is.null(targets_union)) vals <- setdiff(vals, targets_union)
      vals
    }) %>% set_names(sources_of_interest)
  } else {
    result <- list_targets_use
  }
  
  return(result)
}
