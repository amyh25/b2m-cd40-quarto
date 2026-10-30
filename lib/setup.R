# setup.R

library(here)

# paths --------

data_dir <- here("data")
seq1_dir <- file.path(data_dir, "cellranger_output_seq1")
seq2_dir <- file.path(data_dir, "cellranger_output_seq2")

results_dir <- here("results")

# packages ----------

# essential packages
library(Matrix)
library(Seurat)
library(tidyverse)

# io packages
library(readxl)

# plotting packages
library(ggpubr)
library(scattermore)
library(ggridges)
library(cowplot)
library(viridis)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# analysis packages
library(cluster)
library(clustree)

# other
library(glue)
library(ggtext)
library(see)

# custom packages
library(RFunctions)
library(Rsc)
