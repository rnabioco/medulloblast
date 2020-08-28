library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(readxl)
library(here)
library(scbp) # install_github("rnabioco/scbp")
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


get_marker_summaries <- function(so,
                             group_col = "seurat_clusters",
                             outdir = "./",
                             prefix = NULL,
                             min_pct = 0.10,
                             pvalue = 0.05,
                             only_pos = TRUE,
                             tsv_output_dir = NULL,
                             xlsx_output_dir = NULL){

  full_markers <- presto::wilcoxauc(so, group_col)

  mkrs <- filter(full_markers, padj < pvalue,  pct_in > min_pct)
  if(only_pos){
    mkrs <- filter(mkrs, logFC > 0)
  }
  mkrs <- mkrs %>%
    group_by(group) %>%
    arrange(padj,
            desc(logFC),
            .by_group = TRUE)

  if(is.null(prefix)){
    prefix <- ""
  } else {
    prefix <- str_c(prefix, "_")
  }

  if(is.null(tsv_output_dir)){
    tsv_output_dir <- outdir
  }

  if(is.null(xlsx_output_dir)){
    xlsx_output_dir <- outdir
  }

  mkrs %>%
    write_tsv(file.path(tsv_output_dir,
                        str_c(prefix, "cluster_markers.tsv")))

  mkrs %>%
    ungroup() %>%
    split(., .$group) %>%
    write_markers_xlsx(.,
                       file.path(xlsx_output_dir,
                                 str_c(prefix, "cluster_markers.xlsx")))

  dplyr::lst(full_markers, mkrs)
}



get_alra_assay <- function(so, file_name, overwrite = FALSE){

  ## only used to add to cellbrowser
  if(overwrite || !file.exists(file_name)){
    so <- RunALRA(so, setDefaultAssay = FALSE)
    gc()
    alra_assay <- so@assays$alra
    qs::qsave(alra_assay, file_name)
  } else {
    alra_assay <- qs::qread(file_name)
  }
  alra_assay
}
