# gg_seu <- AddOrthoAssay(gg_seu, orthos_comb, assay_type = "RNA", origName = "gene", newName = "zebraFinch_Symbol")
library(readr)
chicken_mouse_zebraFinch_human_orthologs_20240620 <- read_csv(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/geneLists/chicken_mouse_zebraFinch_human_orthologs_20240620.csv")

AddOrthoAssay <- function(seuObj, ortho_table, assay_type, origName, newName) {
  # Create new assay for orthos
  seuObj_rename <- GetAssayData(seuObj, assay = assay_type)
  original_rownames <- rownames(seuObj_rename)
  seuObj_rename <- seuObj_rename[grepl("ENS", rownames(seuObj_rename)), ]
  
  # Change rownames
  ensembl_ids <- ortho_table[[origName]]
  symbols <- ortho_table[[newName]]
  
  matching_symbols <- symbols[match(rownames(seuObj_rename), ensembl_ids)]
  new_rownames <- ifelse(is.na(matching_symbols), original_rownames, matching_symbols) # If there is no match, keep the original name
  
  rownames(seuObj_rename) <- new_rownames
  
  # Add new assay
  new_assay <- CreateAssayObject(counts = seuObj_rename)
  new_assay_name <- paste("Orthos", assay_type, sep = "_")
  seuObj[[new_assay_name]] <- new_assay
  
  # Also keep variable features
  variable_features <- VariableFeatures(seuObj, assay = assay_type)
  
  # Match variable features to symbols
  matching_symbols_vf <- symbols[match(variable_features, ensembl_ids)]
  new_variable_features <- ifelse(is.na(matching_symbols_vf), variable_features, matching_symbols_vf)
  
  # Filter to keep only those that are in the new assay
  new_variable_features <- new_variable_features[new_variable_features %in% rownames(seuObj_rename)]
  new_variable_features <- new_variable_features[!grepl("ENSGALG", new_variable_features)]
  
  VariableFeatures(seuObj[[new_assay_name]]) <- new_variable_features
  
  return(seuObj)
}
