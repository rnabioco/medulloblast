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

#' write out top n markers
#'
top_marker_matrix <-  function(mkrs,
                               n = NULL,
                               min_pct = 0.10,
                               pvalue = 0.05,
                               only_pos = TRUE,
                               col_to_keep = "logFC"){
  to_keep <- mkrs %>%
    filter(padj < pvalue,
           pct_in > min_pct) %>%
    group_by(group) %>%
    arrange(padj, logFC, .by_group = TRUE)

  if(!is.null(n)){
    to_keep <- to_keep %>%
      slice(1:n)
  }

  if(only_pos){
    to_keep <- filter(to_keep, logFC > 0)
  }

  features_to_keep <- to_keep %>%
    pull(feature) %>%
    unique()

  res <- filter(mkrs, feature %in% features_to_keep) %>%
    select(feature, group, !!sym(col_to_keep)) %>%
    pivot_wider(names_from = "group",
              values_from = all_of(col_to_keep))

  res
}
