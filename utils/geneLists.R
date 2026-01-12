library(readxl)

CanMarkers <- 
  c("SNAP25", # neurons
    "SLC17A6", #glut cells,
    # "SATB2", 
    # "SLC1A3",
    "SOX4", 
    "DCX", #migrating neuroblasts
    "NECTIN3",
    "GAD1", "GAD2", #gaba cells
    "SPEF2", # epen
    "NR2E1","PAX6", #RGC
    "LRIG1", "SLC1A2", #Astro
    "PDGFRA", "VCAN", #OPCs
    # "FYN", # COP
    "TCF7L2", #NFOL
    "MBP", "PLP1", #Oligos
    "CSF1R", "C1QB", #Microglia
    "FLT1", "EGFL7" #Endothelial
  )

# Transcriptions genes based on 
TFs <- readxl::read_excel("/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/geneLists/tf_list.xlsx",
                          sheet = "Table S1. Related to Figure 1B", skip = 1) %>% filter(`[[`(., 4) == "Yes") %>% pull(Name)

# Allen iterative gene used
iterAllnGenes <- readRDS(file =  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_AnalysisR/06_Taxonomy/AllenClustering/RDS_files/markersIterClust.rds")



