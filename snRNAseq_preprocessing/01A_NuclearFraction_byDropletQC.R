# Install packages
# devtools::install_github("powellgenomicslab/DropletQC")

#---

# Load packages
library(DropletQC)
library(tidyverse)

# SET WORK FOLDER AS THE WORKING DIRECTORY
setwd("/project/Neuroinformatics_Core/Roberts_lab/s433904/TelEvo/Published_HVCandRA/RAW_for_CUS32_7355_PubHVCRA/01-analysis/03-count")
root_dir <- "/project/Neuroinformatics_Core/Roberts_lab/s433904/TelEvo/Published_HVCandRA/RAW_for_CUS32_7355_PubHVCRA/01-analysis/03-count"

# Input list
InputList <- basename(list.dirs(recursive = FALSE))

for(x in 1:length(InputList)){
  nf1 <- 
    nuclear_fraction_tags(
      outs = file.path(root_dir, sprintf("%s/outs", InputList[x])),
      cores = 8,
      verbose = T
    )

  save(nf1, file = file.path(root_dir, sprintf("%s/outs/%s_IntronicReadRatios.Rda", InputList[x], InputList[x])))

}


