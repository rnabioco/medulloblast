library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(tidyverse)
library(readxl)
library(here)
library(scbp)
library(presto)
library(qs)
theme_set(theme_cowplot())


# ----------------------------------------------------
library(future)

plan("multicore", workers = 7)
options(future.globals.maxSize = 2 * 1024 ^ 3)

proj_dir <- here()
data_dir <- file.path(proj_dir, "data", "cellranger", "results")
doc_dir <- file.path(proj_dir, "docs")

fig_dir <- "figs"
mkrs_dir <- "markers"
tbls_dir <- "tables"
obj_dir <- "objects"
walk(c(fig_dir, mkrs_dir, tbls_dir, obj_dir),
     dir.create, showWarnings = F)

out_dirs <- dplyr::lst(fig_dir, mkrs_dir, tbls_dir, obj_dir)
