library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(tidyverse)
library(readxl)
library(here)
library(scbp)
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
xcel_dir <- file.path(mkrs_dir, "xlsx")
walk(c(fig_dir, mkrs_dir, tbls_dir, xcel_dir),
     dir.create,
     showWarnings = F)
