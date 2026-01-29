# LOAD PACKAGES
library(eulerr)

# DIRECTORIES
Fig_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/21_DGE_whyNeuronsSpecial/Figures"

# LOAD SOME GENES
# ASD genes
sfari_asd <- read_csv("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/geneLists/SFARI-Gene_genes_04-03-2025release_04-22-2025export.csv")
asd_genes <- sfari_asd$`gene-symbol` %>% unique()

# By TFs
tf_list <- readxl::read_excel("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/geneLists/tf_list.xlsx",
  sheet = "Table S1. Related to Figure 1B",
  skip = 1
)

colnames(tf_list)[colnames(tf_list) == "...4"] <- "isTF"
TFs <- tf_list$Name[tf_list$isTF == "Yes"]


# VISUALIZE W/SINGLE CELL
## set up
gg_gabaOfI_comb <- RenameAssays(gg_gabaOfI, SCT_ortho = "SCTcmb")
seuObj_gaba_comb <- RenameAssays(seuObj_gaba, SCTorigIdent3 = "SCTcmb")

common_features <- intersect(
  rownames(gg_gabaOfI_comb[["SCTcmb"]]),
  rownames(seuObj_gaba_comb[["SCTcmb"]])
)
gg_gabaOfI_comb <- subset(gg_gabaOfI_comb, features = common_features)
seuObj_gaba_comb <- subset(seuObj_gaba_comb, features = common_features)

gg_gabaOfI_comb@assays[["RNA_ortho"]] <- NULL
seuObj_gaba_comb@assays[["RNA"]] <- NULL

combGABA <- merge(gg_gabaOfI_comb, seuObj_gaba_comb)

combGABA$clusterComb <- ifelse(is.na(combGABA$cn3_subclass),
  "gg_PV_SST", combGABA$cn3_subclass
)

## visualize
### venn diagram
fit <- euler(list(
  SSTsong = sst_common_up_spec,
  PVsong = pv_common_up_spec
))

# Plot with numbers inside
venDiag <- plot(fit,
  quantities = TRUE, # show counts
  fills = c("skyblue", "salmon"),
  labels = list(col = "black", font = 2)
)

ggsave(
  filename = file.path(Fig_dir, paste0("venndiagr_", "song-mge", ".svg")), plot = venDiag,
  width = 4, height = 2
)


### pick interesting genes
asd_tfs <- intersect(TFs, asd_genes)

pv_song_GOI_asdTF <- setdiff(intersect(pv_common_up_spec, asd_tfs), up_both) # this is to keep in mind
sst_song_GOI_asd_TF <- setdiff(intersect(sst_common_up_spec, asd_tfs), up_both)

pv_song_GOI_asd <- setdiff(intersect(pv_common_up_spec, asd_genes), up_both) # this is to keep in mind
sst_song_GOI_asd <- setdiff(intersect(sst_common_up_spec, asd_genes), up_both)

pv_song_GOI_TF <- setdiff(intersect(pv_common_up_spec, TFs), up_both) # this is to keep in mind
sst_song_GOI_TF <- setdiff(intersect(sst_common_up_spec, TFs), up_both)

pv_song_GOI <- setdiff(intersect(pv_common_up_spec, c(TFs, asd_genes)), up_both)
sst_song_GOI <- setdiff(intersect(sst_common_up_spec, c(TFs, asd_genes)), up_both)


### make plots
geneLists <- list(
  up_both = up_both,
  pv_song_picked = pv_song_GOI, # pv_common_up_spec or pv_song_GOI
  sst_song_picked = sst_song_GOI # sst_common_up_spec or sst_song_GOI
)

for (nm in names(geneLists)) {
  use <- geneLists[[nm]]

  p <- Stacked_VlnPlot(
    combGABA,
    features = use, group.by = "clusterComb", x_lab_rotate = TRUE
  ) & scale_x_discrete(limits = c("gg_PV_SST", "MGE_song", "SST", "PVALB", "PVALB_song")) &
    theme(
      axis.line = element_line(color = "#060606")
    )

  ggsave(
    filename = file.path(Fig_dir, paste0("vlnplot_", nm, ".svg")), plot = p,
    width = 4, height = length(use) / 2
  )
}

