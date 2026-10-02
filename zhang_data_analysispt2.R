# load packages and set up environment
library(Seurat)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(data.table)
library(stringr)
setwd("/gpfs1/pi/kmccreer/scRNAseq/issy_26-27/honors_thesis")

## Visualizing data! 
obj <- readRDS("merged_seurat_obj_allRNAsamples.rds"")