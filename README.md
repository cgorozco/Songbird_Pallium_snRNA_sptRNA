# Songbird Pallium single-nuclei RNA and spatial transcriptomics
Cellular and spatial transcriptomics of several regions across the zebra finch pallium; Hyperpallium, CM, AM, HVC, NCM, NCL, LMAN, AN, RA, and AA. Providing a transcriptomic atlas of both evolutionarily novel and conserved brain regions of the bird pallium (dorsal telencephalon), spanning areas analogous to the mammalian motor, premotor, prefrontal, secondary auditory, and primary visual cortices. We focused on identifying potential cellular innovations in song-learning brain regions that underlie the eponymous vocal imitation behavior, as these regions and behavior represent evolutionarily novel specializations in the avian lineage.

## Data
Seurat processed data objects can be accessed at https://cloud.biohpc.swmed.edu/______
Raw and processed data are available in the Gene Expression Omnibus for the snRNAseq (GSE316328) and for the spatial dataset (GSE316807). 
To run the analysis, directories will need to be changed.

## Directories
comparative - Comparisons between songbird and chicken neurons. Scripts are divided by excitatory or inhibitory neurons at different resolutions of clustering

snRNAseq_analysis – Subdirectories for downstream analyses of the processed zebra finch pallium snRNA-seq dataset. Including LISI scoring, classifying excitatory neurons of song-sorround brain regions, LIANA-based cell–cell communication, ...

snRNAseq_preprocessing - Processing of snRNA-seq data following CellBender (ambient RNA removal). Includes further processing for nuclear fraction, doublets, batch effects, merging, integration, and itertive consensus clustering (scrattch.hicat).

spatial_RNAseq - Analysis of spatial transcriptomic dataset for clustering, processing, and integrating with snRNAseq dataset

paper figures - Script for making portions of the figures in paper using the zebra finch snRNAseq dataset

utils - Several files with utility functions used in other scripts
