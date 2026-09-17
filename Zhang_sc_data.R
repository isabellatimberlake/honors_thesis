### Data exploration of scRNAseq from Zhang et al. for thesis

# install.packages("Seurat")
# load packages
library(Seurat)
library(tidyverse)
setwd("/gpfs1/pi/kmccreer/scRNAseq/issy_26-27/issy_thesis")

# load data
PT653 <- Read10X("/gpfs1/pi/kmccreer/scRNAseq/issy_26-27/issy_thesis/GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz")
head(PT653)

