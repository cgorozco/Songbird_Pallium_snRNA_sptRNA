library(readxl)

seuObj$clusterName3 <- factor(seuObj$clusterName3, levels = dendro_order_3)
Idents(seuObj) <- "clusterName3"

# add coarse lvl
coarse_lvl <- readxl::read_excel(
  "/project/Neuroinformatics_Core/Roberts_lab/s433904/PalliumEvo_paper/Tables/CellClassesClassification.xlsx", 
  sheet = "ForCode")

metaUse <- FetchData(seuObj, vars = c("clusterName3"))
metaUse$cn3_subclass <- coarse_lvl$Subclass[
  match(metaUse$clusterName3, coarse_lvl$Cluster)]
metaUse$cn3_class <- coarse_lvl$Class[
  match(metaUse$clusterName3, coarse_lvl$Cluster)]

metaOrder <- metaUse %>%
  distinct() %>% 
  arrange(clusterName3, dendro_order_3) %>% 
  pull(cn3_subclass) %>% 
  unique()

metaUse$cn3_subclass <- factor(metaUse$cn3_subclass, levels = metaOrder)

seuObj <- AddMetaData(seuObj, metadata = metaUse)
