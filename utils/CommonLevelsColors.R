class_lvl1_order <- c("Glutamatergic", "GABAergic", "Neurogenic",
                      "Astro", "Epen", "OPC_Oligo", "Micro", "Endo")

#-- orig.ident levels
custom_colors <- c(
  colorRampPalette(c("#663399", "#d8b3ff"))(4),
  colorRampPalette(c("#336666", "#99cccc"))(3),
  colorRampPalette(c("#336600", "#99cc66"))(3),
  colorRampPalette(c("#cc9900", "#ffff00"))(4),
  colorRampPalette(c("#000099", "#6699ff"))(3),
  colorRampPalette(c("#333333", "#cccccc"))(5),
  colorRampPalette(c("#660000", "#ff0000"))(4),
  colorRampPalette(c("#006666", "#00cccc"))(3),
  colorRampPalette(c("#663333", "#ffcccc"))(5),
  colorRampPalette(c("#cc3000", "#ff9966"))(4)
)

origIdent_colors <- c(
  colorRampPalette(c("#333333", "#cccccc"))(4), # hyper
  colorRampPalette(c("#553821", "#A1846E"))(3),
  colorRampPalette(c("#4E472F", "#9B937B"))(3),
  colorRampPalette(c("#163E29", "#628B76"))(4),
  colorRampPalette(c("#3B5637", "#87A384"))(3),
  colorRampPalette(c("#4E6A76", "#9AB6C2"))(5),
  colorRampPalette(c("#32475F", "#7F94AC"))(4),
  colorRampPalette(c("#323257", "#7E7EA3"))(3),
  colorRampPalette(c("#8C3A29", "#D88775"))(5),
  colorRampPalette(c("#93844E", "#E0D09A"))(4)
)

region_colors <- c("#6F6F6F", "#7A5030", "#706543", "#1F593B", "#547B4F", "#6F97A8", "#486688", "#47477C", "#C8533A", "#D2BC6F")

# items <- unique(stackMeta$orig.ident)
# sorted_items <- items[order(factor(gsub("_.*", "", items), levels = lvl_order))]
origIdent_lvls <- c("Hyper_1", "Hyper_2", "Hyper_3", "Hyper_DPM", "Meso_2", "Meso_3", "Meso_DPM",
                    "Av_1", "Av_2", "Av_3", "LMAN_1", "LMAN_2", "LMAN_3", "LMAN_DPM", 
                    "AN_1", "AN_2", "AN_DPM", "HVC_2", "HVC_3", "HVC_DPM", "HVC_Pub1", "HVC_Pub2",
                    "NC_1", "NC_2", "NC_3", "NC_DPM", "NCM_1", "NCM_3", "NCM_DPM",
                    "RA_1", "RA_2", "RA_3", "RA_DPM", "RA_Pub1", "Arco_1", "Arco_2", "Arco_3", "Arco_DPM")


#-- colors by region
lvl_order <- c("Hyper", "Meso", "Av", "LMAN", "An", "HVC", "NC", "NCM", "RA", "Arco")
get_middle_color <- function(palette_func, n_colors) {
  colors <- palette_func(n_colors)
  middle_index <- ceiling(n_colors / 2)  # Get the middle index
  return(colors[middle_index])
}
custom_colors_rgn <- c(
  get_middle_color(colorRampPalette(c("#663399", "#d8b3ff")), 5),  # For "Hyper"
  get_middle_color(colorRampPalette(c("#336666", "#99cccc")), 3),  # For "Meso"
  get_middle_color(colorRampPalette(c("#336600", "#99cc66")), 3),  # For "Av"
  get_middle_color(colorRampPalette(c("#cc9900", "#ffff00")), 4),  # For "LMAN"
  get_middle_color(colorRampPalette(c("#000099", "#6699ff")), 3),  # For "An"
  get_middle_color(colorRampPalette(c("#333333", "#cccccc")), 5),  # For "HVC"
  get_middle_color(colorRampPalette(c("#660000", "#ff0000")), 4),  # For "NC"
  get_middle_color(colorRampPalette(c("#006666", "#00cccc")), 3),  # For "NCM"
  get_middle_color(colorRampPalette(c("#663333", "#ffcccc")), 5),  # For "RA"
  get_middle_color(colorRampPalette(c("#cc3000", "#ff9966")), 4)   # For "Arco"
)
color_mapping <- setNames(custom_colors_rgn, lvl_order)
PallRgn_cols <- c("#9b59b6", "#669966", "#6699cc", "#e38d74")

custom_colors_rgn2 <- c("#6F6F6F", "#7A5030FF", "#706543FF", "#1F593BFF", "#547B4FFF", "#6F97A8FF", "#486688FF", "#47477CFF", "#C8533AFF", "#D2BC6FFF")

#-- levels of chicken glut levels
gg.levels.lvl1 <- c("Ex.DACH2.CALCR", "Ex.DACH2.SV2C", "Ex.SATB2", "Ex.SATB2.KIAA1217", "Ex.CACNA1H",  "Ex.TCF7L2")
gg_levels_lvl1 <- c("Ex_DACH2_CALCR", "Ex_DACH2_SV2C", "Ex_SATB2", "Ex_SATB2_KIAA1217", "Ex_CACNA1H",  "Ex_TCF7L2")
gg_levels_lvl1_v2 <- c("Ex_CACNA1H", "Ex_DACH2_CALCR", "Ex_DACH2_SV2C", "Ex_SATB2", "Ex_SATB2_KIAA1217",  "Ex_TCF7L2")
gg.levels.lvl1.v2 <- c("Ex.CACNA1H", "Ex.DACH2.CALCR", "Ex.DACH2.SV2C", "Ex.SATB2", "Ex.SATB2.KIAA1217",  "Ex.TCF7L2")

gg_levels_lvl3 <- c("Ex_TSHZ2_NR4A2", # hyper
                    "Ex_DACH2_ITGA9", "Ex_DACH2_GRIK4", "Ex_DACH2_ZMAT4", "Ex_DACH2_TAC1",  # hyper/nido
                    "Ex_DACH2_ADAMTS5", "Ex_DACH2_MGAT4C", "Ex_DACH2_NR4A3", "Ex_DACH2_LHX2",
                    "Ex_DACH2_CEMIP", "Ex_DACH2_LUZP2", "Ex_DACH2_RORB", "Ex_DACH2_SLIT2",
                    "Ex_Pre_SATB2", # Meso other
                    "Ex_SATB2_FOXP2", "Ex_SATB2_ZNF385B", "Ex_SATB2_SOX6", "Ex_SATB2_OVOA", 
                    "Ex_KIAA1217", "Ex_KIAA1217_BCL6", #Meso other
                    "Ex_CACNA1H_LHX9", "Ex_CACNA1H_PROX1", "Ex_CACNA1H_CPA6", "Ex_CACNA1H_KIT", "Ex_CACNA1H_MCTP2", # arco and DG
                    "Ex_TCF7L2", # cluster w/gaba
                    "Ex_BCL6", # peri and/or olf
                    "Ex_Pre_KCNH7" # other pre
                    )

gg.levels.lvl3 <- c("Ex.TSHZ2.NR4A2", # hyper
                    "Ex.DACH2.ITGA9", "Ex.DACH2.GRIK4", "Ex.DACH2.ZMAT4", "Ex.DACH2.TAC1",  # hyper/nido
                    "Ex.DACH2.ADAMTS5", "Ex.DACH2.MGAT4C", "Ex.DACH2.NR4A3", "Ex.DACH2.LHX2",
                    "Ex.DACH2.CEMIP", "Ex.DACH2.LUZP2", "Ex.DACH2.RORB", "Ex.DACH2.SLIT2",
                    "Ex.Pre.SATB2", # Meso other
                    "Ex.SATB2.FOXP2", "Ex.SATB2.ZNF385B", "Ex.SATB2.SOX6", "Ex.SATB2.OVOA", 
                    "Ex.KIAA1217", "Ex.KIAA1217.BCL6", #Meso other
                    "Ex.CACNA1H.LHX9", "Ex.CACNA1H.PROX1", "Ex.CACNA1H.CPA6", "Ex.CACNA1H.KIT", "Ex.CACNA1H.MCTP2", # arco and DG
                    "Ex.TCF7L2", # cluster w/gaba
                    "Ex.BCL6", # peri and/or olf
                    "Ex.Pre.KCNH7" # other pre
                    )

gg_levels_lvl3_v2 <- c("Ex_CACNA1H_LHX9", "Ex_CACNA1H_PROX1", "Ex_CACNA1H_CPA6", "Ex_CACNA1H_KIT", "Ex_CACNA1H_MCTP2", # arco and DG
                       "Ex_DACH2_ITGA9", "Ex_DACH2_GRIK4", "Ex_DACH2_ZMAT4", "Ex_DACH2_TAC1",  # hyper/nido
                       "Ex_DACH2_ADAMTS5", "Ex_DACH2_MGAT4C", "Ex_DACH2_NR4A3", "Ex_DACH2_LHX2",
                       "Ex_DACH2_CEMIP", "Ex_DACH2_LUZP2", "Ex_DACH2_RORB", "Ex_DACH2_SLIT2",
                       "Ex_TCF7L2", # cluster w/gaba
                       "Ex_TSHZ2_NR4A2", # hyper
                       "Ex_KIAA1217", "Ex_KIAA1217_BCL6", #Meso other
                       "Ex_SATB2_ZNF385B", "Ex_SATB2_SOX6", "Ex_SATB2_FOXP2", "Ex_SATB2_OVOA",
                       "Ex_BCL6" # peri and/or olf
                       )

gg.levels.lvl3.v2 <- gsub("_", ".", gg_levels_lvl3_v2)

###==== new as of 1/17/2025
dendro_order_3 <- c("Glut-LMAN(44)", "Glut-LMAN(45)", "Neurogenic-NidoHyperPall(42)", "Neurogenic(4)", 
                    "Glut-RA(19)", "Glut-NidoHyperPall(23)", "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)",
                    "Glut-MesoPall(6)", "Neurogenic-NidoHyperPall(41)", "Neurogenic-NidoHyperPall(43)",
                    "Glut-MesoPall(5)", "Glut-MesoPall(29)", "Glut-MesoPall(30)", "Glut-Hyper(9)",
                    "Glut-Arco(1)", "Glut-RA(18)", "Glut-Arco(10)", "Glut-ArcoPall(27)", "Glut-MesoPall(2)",
                    "Glut(21)", "Glut-NidoHyperPall(20)", "Glut-NidoHyperPall(25)", "Glut-NidoHyperPall(24)",
                    "Glut-NidoHyperPall(26)", "LGE_Pre(35)", "LGE(3)", "MGE_ST18(52)", "MGE_ST18(50)", "MGE_ST18(51)",
                    "MGE_PVALB(40)", "MGE_PVALB(31)", "MGE_PVALB(34)", "MGE_SST(32)", "MGE_SST(36)", "MGE_SST(37)",
                    "MGE_SST(38)", "MGE_song(33)", "MGE_song(39)", "CGE(11)", "CGE_LHX8(47)", "CGE_LHX8(48)",
                    "CGE_LHX8(49)", "CGE_Lamp5(13)", "Epen-NidoHyperPall(17)", "Epen-NidoHyperPall(14)", 
                    "Astro(12)", "OPC(53)", "COP-NFOL(58)", "COP-NFOL(55)", "COP-NFOL(54)", "Oligo(57)", 
                    "Oligo(56)", "Micro(16)", "Endo(15)") # with chatGPT help

dendro.order.3 <- gsub("[_()-]", ".", dendro_order_3)


glutLvls_byRgnPct <- c("Glut-RA(19)", "Glut-RA(18)", "Glut-ArcoPall(27)", "Glut-Arco(1)", "Glut-Arco(10)",
                       "Glut-HVC(7)", "Glut-HVC(8)", "Glut-NidoHyperPall(28)", "Glut-NidoHyperPall(23)",
                       "Glut-NidoHyperPall(24)", "Glut-NidoHyperPall(26)", "Glut-NidoHyperPall(20)",
                       "Glut-LMAN(45)", "Glut-LMAN(44)", 
                       "Glut(21)", 
                       "Glut-NidoHyperPall(25)", "Glut-Hyper(9)",
                       "Glut-MesoPall(6)", "Glut-MesoPall(5)", "Glut-MesoPall(29)", "Glut-MesoPall(30)",  "Glut-MesoPall(2)")

glutLvls.byRgnPct <- gsub("[()\\-]", ".", glutLvls_byRgnPct)

##---- Human rois
# human_roiLvls <- c("A13", "A14", "A19", "A1C", "A23", "A24", "A25", "A29-A30", "A32", "A35-A36", 
#   "A35r", "A38", "A40", "A43", "A44-A45", "A46", "A5-A7", "FI", "Idg", 
#   "Ig", "ITG", "LEC", "M1C", "MEC", "MTG", "Pro", "S1C", "STG", "TF", 
#   "TH-TL", "V1C", "V2", 
#   "CA1C-CA3C", "CA1R-CA2R", "CA1R-CA2R-CA3R", "CA1U", "CA1U-CA2U-CA3U", # hippocampus
#   "CA2U-CA3U", "CA3R", "CA4C-DGC", "DGR-CA4Rpy", "DGU-CA4Upy", "Sub", 
#   "Cla", # lateral pall
#   "Pir", "AON", # olfactory regs
#   "BL", "BM", "CEN", "CMN", "CoA", "La", # amygdala
#   
#   "BNST",
#   "SEP", "SI", "CaB", "GPe", "GPi", "NAC", "Pu",  
#   "CbDN", "CBL", "CBV",
#   "HTHma", "HTHma-HTHtub", "HTHpo", "HTHpo-HTHso", "HTHso", "HTHso-HTHtub", "HTHtub", "MN",
#   "ANC", "CM", "CM-Pf", "ETH", "LG", "LP", "LP-VPL", "MD", "MD-Re", "MG", "Pul", "STH", "VA", "VLN", "VPL",
#   
#   "IC", "PAG", "PAG-DR", "PTR", "RN", "SC", "SN", "SN-RN", 
#   "IO", "MoAN", "MoRF-MoEN", "MoSR", 
#   "DTg", "PB", "PN", "PnAN", "PnEN", "PnRF",
#   "SpC"
# )











#####========= ARCHIVE ====================================================#####

# new_levels_dot <- c("Glutamatergic.18", "Glutamatergic.21", "Glutamatergic.11", "Glutamatergic.20", "Glutamatergic.30",
#                     "Glutamatergic.16", "Glutamatergic.3", "Glutamatergic.2", "Glutamatergic.0", "Glutamatergic.7",
#                     "Glutamatergic.4", "Glutamatergic.22", "Glutamatergic.32", "Glutamatergic.26", "Neurogenic.13", 
#                     "Neurogenic.15", "Neurogenic.34", "Neurogenic.37", "LGE.14", "LGE.10", "MGE.Sst.31", "MGE.Pvalb.8",
#                     "MGE.Sst.17", "MGE.Sst.24", "MGE.Sst.36", "CGE.Lamp5.27", "CGE.Vip.25", "CGE.Vip.Sncg.29", "CGE.Lamp5.33",
#                     "Glutamatergic.35", "Epen.23", "Astro.1", "Astro.9", "Astro.6", "OPC.12", "TOL.38", "Oligo.5", "Micro.28",
#                     "Endo.19")
# 
# new_levels <- c("Glutamatergic_18", "Glutamatergic_21", "Glutamatergic_11", "Glutamatergic_20", "Glutamatergic_30",
#                 "Glutamatergic_16", "Glutamatergic_3", "Glutamatergic_2", "Glutamatergic_0", "Glutamatergic_7",
#                 "Glutamatergic_4", "Glutamatergic_22", "Glutamatergic_32", "Glutamatergic_26", "Neurogenic_13", 
#                 "Neurogenic_15", "Neurogenic_34", "Neurogenic_37", "LGE_14", "LGE_10", "MGE_Sst_31", "MGE_Pvalb_8",
#                 "MGE_Sst_17", "MGE_Sst_24", "MGE_Sst_36", "CGE_Lamp5_27", "CGE_Vip_25", "CGE_Vip_Sncg_29", "CGE_Lamp5_33",
#                 "Glutamatergic_35", "Epen_23", "Astro_1", "Astro_9", "Astro_6", "OPC_12", "TOL_38", "Oligo_5", "Micro_28",
#                 "Endo_19")

# #-- renaming of gluts and ordering 
# # Create a named vector for the mappings
# glut_rename <- c(
#   "Glutamatergic_22" = "Pan_glut(22)",
#   "Glutamatergic_0" = "NonArco_glut(0)",
#   "Glutamatergic_2" = "NonArco_glut(2)",
#   "Glutamatergic_16" = "NonArco_glut(16)",
#   "Glutamatergic_3" = "NonArco_glut(3)",
#   "Glutamatergic_35" = "LMAN_glut(35)",
#   "Glutamatergic_11" = "HVC_glut(11)",
#   "Glutamatergic_21" = "Meso_glut(21)",
#   "Glutamatergic_7"  = "Meso_glut(7)",
#   "Glutamatergic_20" = "Meso_glut(20)",
#   "Glutamatergic_32" = "Meso_glut(32)",
#   "Glutamatergic_26" = "Meso_glut(26)",
#   "Glutamatergic_30" = "Meso_glut(30)",
#   "Glutamatergic_4"  = "Arco_glut(4)",
#   "Glutamatergic_18" = "RA_glut(18)"
# )
# 
# # Define the levels in the desired order
# glut_levels <- c(
#   "Pan_glut(22)", "NonArco_glut(0)", "NonArco_glut(2)", "NonArco_glut(16)", "NonArco_glut(3)", 
#   "LMAN_glut(35)", "HVC_glut(11)", "Meso_glut(21)", "Meso_glut(7)", "Meso_glut(20)", 
#   "Meso_glut(32)", "Meso_glut(26)", "Meso_glut(30)", "Arco_glut(4)", "RA_glut(18)"
# )
# 
# # # Apply the renaming and specify the levels for a column
# # df$your_column <- factor(recode(df$your_column, !!!glut_rename), levels = glut_levels)
# # rownames(df) <- factor(recode(rownames(df), !!!glut_rename), levels = glut_levels)
# 
# glut.rename <- c(
#   "Glutamatergic.22" = "Pan.glut(22)",
#   "Glutamatergic.0" = "NonArco.glut(0)",
#   "Glutamatergic.2" = "NonArco.glut(2)",
#   "Glutamatergic.16" = "NonArco.glut(16)",
#   "Glutamatergic.3" = "NonArco.glut(3)",
#   "Glutamatergic.35" = "LMAN.glut(35)",
#   "Glutamatergic.11" = "HVC.glut(11)",
#   "Glutamatergic.21" = "Meso.glut(21)",
#   "Glutamatergic.7"  = "Meso.glut(7)",
#   "Glutamatergic.20" = "Meso.glut(20)",
#   "Glutamatergic.32" = "Meso.glut(32)",
#   "Glutamatergic.26" = "Meso.glut(26)",
#   "Glutamatergic.30" = "Meso.glut(30)",
#   "Glutamatergic.4"  = "Arco.glut(4)",
#   "Glutamatergic.18" = "RA.glut(18)"
# )
# 
# # Define the levels in the desired order
# glut.levels <- c(
#   "Pan.glut(22)", "NonArco.glut(0)", "NonArco.glut(2)", "NonArco.glut(16)", "NonArco.glut(3)", 
#   "LMAN.glut(35)", "HVC.glut(11)", "Meso.glut(21)", "Meso.glut(7)", "Meso.glut(20)", 
#   "Meso.glut(32)", "Meso.glut(26)", "Meso.glut(30)", "Arco.glut(4)", "RA.glut(18)"
# )
# 
# # same as above but for all cells
# all_rename<- c(
#   "Glutamatergic_22" = "Pan_glut(22)",
#   "Glutamatergic_0" = "NonArco_glut(0)",
#   "Glutamatergic_2" = "NonArco_glut(2)",
#   "Glutamatergic_16" = "NonArco_glut(16)",
#   "Glutamatergic_3" = "NonArco_glut(3)",
#   "Glutamatergic_35" = "LMAN_glut(35)",
#   "Glutamatergic_11" = "HVC_glut(11)",
#   "Glutamatergic_21" = "Meso_glut(21)",
#   "Glutamatergic_7"  = "Meso_glut(7)",
#   "Glutamatergic_20" = "Meso_glut(20)",
#   "Glutamatergic_32" = "Meso_glut(32)",
#   "Glutamatergic_26" = "Meso_glut(26)",
#   "Glutamatergic_30" = "Meso_glut(30)",
#   "Glutamatergic_4"  = "Arco_glut(4)",
#   "Glutamatergic_18" = "RA_glut(18)",
#   "Neurogenic_13" = "NGC_glut(13)",
#   "Neurogenic_15" = "NGC_glut(15)",
#   "Neurogenic_34" = "NGC_gaba(34)",
#   "Neurogenic_37" = "NGC_gaba(37)"
# )
# 
# all_levels <- c(
#   "Pan_glut(22)", "NonArco_glut(0)", "NonArco_glut(2)", "NonArco_glut(16)", "NonArco_glut(3)", 
#   "LMAN_glut(35)", "HVC_glut(11)", "Meso_glut(21)", "Meso_glut(7)", "Meso_glut(20)", 
#   "Meso_glut(32)", "Meso_glut(26)", "Meso_glut(30)", "Arco_glut(4)", "RA_glut(18)",
#   "NGC_glut(13)", "NGC_glut(15)", "NGC_gaba(34)", "NGC_gaba(37)", 
#   "LGE_14", "LGE_10", "MGE_Sst_31", "MGE_Pvalb_8",
#   "MGE_Sst_17", "MGE_Sst_24", "MGE_Sst_36", "CGE_Lamp5_27", "CGE_Vip_25", "CGE_Vip_Sncg_29", "CGE_Lamp5_33",
#   "Epen_23", "Astro_1", "Astro_9", "Astro_6", "OPC_12", "TOL_38", "Oligo_5", 
#   "Micro_28",  "Endo_19")


# dendro_order_1 <- c(rev(c("OPC(53)", "Epen-NidoHyperPall(14)", "Epen-NidoHyperPall(17)", "NotPct80(22)",
#                           "Astro(12)", "GABA(11)", "CGE_Lamp5(47)", "CGE_Lamp5(48)", "CGE_Lamp5(49)", "CGE_Lamp5(13)", "LGE(3)", "Neurogenic(35)", "MGE_Sst(52)",
#                           "MGE_Sst(50)", "MGE_Sst(51)", "MGE_Sst(32)", "MGE_Sst-NidoHyperPall(36)", "MGE_Sst(38)",
#                           "MGE_Sst(37)", "MGE_Sst(33)", "MGE_Sst-HVC(39)", "GABA(46)", "MGE_Pvalb(40)",
#                           "MGE_Pvalb(31)", "GABA(34)", "Glut-LMAN(44)",
#                           "Glut-LMAN(45)", "Neurogenic-NidoHyperPall(42)", "Neurogenic(4)", "Glut-RA(19)",
#                           "Glut-NidoHyperPall(23)", "Glut-HVC(7)", "Glut-NidoHyperPall(28)", "Glut-HVC(8)",
#                           "Glut-MesoPall(6)", "Neurogenic-NidoHyperPall(41)", "Neurogenic-NidoHyperPall(43)",
#                           "Glut-MesoPall(5)", "Glut-MesoPall(29)", "Glut-MesoPall(30)", "Glut-Hyper(9)",
#                           "Glut-Arco(1)", "Glut-RA(18)", "Glut-Arco(10)", "Glut-ArcoPall(27)", "Glut-MesoPall(2)",
#                           "Glut(21)", "Glut-NidoHyperPall(20)", "Glut-NidoHyperPall(25)", "Glut-NidoHyperPall(24)",
#                           "Glut-NidoHyperPall(26)")), "TOL(55)", "TOL(54)", "OPC_Oligo(58)", "Oligo(56)",
#                     "Oligo(57)", "Micro(16)", "Endo(15)")
# 
# dendro_order <- c(rev(c("OPC(53)", "Epen-NidoHyperPall(14)", "Epen-NidoHyperPall(17)", "NotPct80(22)", "Astro(12)",
#                         "CGE(11)", "CGE_LHX8(47)", "CGE_LHX8(48)", "CGE_LHX8(49)", "CGE_Lamp5(13)",
#                         "LGE(3)", "LGE_Pre(35)", "MGE_ST18(52)",  "MGE_ST18(50)", "MGE_ST18(51)",
#                         "MGE_SST(32)", "MGE_SST(36)", "MGE_SST(38)", "MGE_SST(37)", "MGE_song(33)", "MGE_song(39)",
#                         "mge_pvalb(46)", "MGE_PVALB(40)", "MGE_PVALB(31)", "MGE_PVALB(34)",
#                         "Glut-LMAN(44)", "Glut-LMAN(45)", "Neurogenic-NidoHyperPall(42)", "Neurogenic(4)", "Glut-RA(19)",
#                         "Glut-NidoHyperPall(23)", "Glut-HVC(7)", "Glut-NidoHyperPall(28)", "Glut-HVC(8)",
#                         "Glut-MesoPall(6)", "Neurogenic-NidoHyperPall(41)", "Neurogenic-NidoHyperPall(43)",
#                         "Glut-MesoPall(5)", "Glut-MesoPall(29)", "Glut-MesoPall(30)", "Glut-Hyper(9)",
#                         "Glut-Arco(1)", "Glut-RA(18)", "Glut-Arco(10)", "Glut-ArcoPall(27)", "Glut-MesoPall(2)",
#                         "Glut(21)", "Glut-NidoHyperPall(20)", "Glut-NidoHyperPall(25)", "Glut-NidoHyperPall(24)",
#                         "Glut-NidoHyperPall(26)")), "TOL(55)", "TOL(54)", "OPC_Oligo(58)", "Oligo(56)",
#                   "Oligo(57)", "Micro(16)", "Endo(15)")

# 
# if (exists("seuObj") && any(colnames(seuObj@meta.data) == "Allen_classLvl2")) {
#   cluster_colors <- custom_colors <- c(
#     sample(colorRampPalette(c("#116834", "#b3c6bb"))(length(unique(seuObj$clusterName[seuObj$Allen_classLvl2 == "Glutamatergic"])))),
#     sample(colorRampPalette(c("#34495e", "#99cccc"))(length(unique(seuObj$clusterName[seuObj$Allen_classLvl2 == "GABAergic"])))),
#     sample(colorRampPalette(c("#ba7b43", "#cdbbac"))(length(unique(seuObj$clusterName[!(seuObj$Allen_classLvl2 %in% c("GABAergic", "Glutamatergic"))]))))
#     )
#   cluster_colors <- setNames(cluster_colors, dendro_order)
# }




