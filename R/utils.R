library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(tidyverse)
library(readxl)
library(here)
theme_set(theme_cowplot())

palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

discrete_palette_default <- c(palette_OkabeIto,
                              brewer.pal(12, "Paired"),
                              brewer.pal(9, "Set1"),
                              brewer.pal(8, "Dark2"))
#' Plot cells in reduced dimensionality 2D space
#'
#' @description Cells can be colored by gene or feature in meta.data dataframe
#'
#' @param seurat_obj object of class Seurat
#' @param feature feature to plot, either gene name or column in seurat_obj@meta.data
#' @param plot_dat supplemental data.frame containing feature to plot.
#' Must have a column named cell that contains matching colnames in seurat_obj@data
#' @param pt_size size of points produced by geom_point
#' @param pt_alpha alpha value for points plotted by geom_point
#' @param label_text if TRUE display feature labels on plot
#' @param label_size size of label text
#' @param label_color color of label text
#' @param .cols vector of colors to use for plot.
#' @param cell_filter character vector of cell names to include in plot
#' @param palette_type color palette type to use (either viridis, brewer, or cloupe)
#' defaults to using cellranger loupe-like colors
#' @param col_pal palette name to use if palette_type is brewer
#' @param max_y maximum feature value to set scale to. Defaults to max of the feature
#' @param legend_title string to supply for title for the legend
#' @param embedding dimensionality reduction to extract from seurat_obj. Can be any
#' dr method present in seurat_obj@dr (e.g. umap, pca, tsne). defaults to tsne
#'
plot_feature <- function(seurat_obj,
                         feature = NULL,
                         plot_dat = NULL,
                         pt_size = 0.001,
                         pt_alpha = 1,
                         label_text = FALSE,
                         label_size = 6,
                         label_color = "grey",
                         .cols = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         legend_title = NULL,
                         embedding = "tsne"){

  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")

  if(!embedding %in% names(seurat_obj@reductions)){
    stop(paste0(embedding, " not found in seurat object"))
  }

  embed_dat <- seurat_obj@reductions[[embedding]]@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  embed_cols <- colnames(embed_dat)
  xcol <- embed_cols[2]
  ycol <- embed_cols[3]

  embed_dat <- left_join(mdata, embed_dat, by = "cell")

  if (!is.null(cell_filter)){
    embed_dat <- dplyr::filter(embed_dat,
                               cell %in% cell_filter)
  }

  meta_data_col <- feature %in% colnames(embed_dat)

  if (!is.null(feature) & !meta_data_col) {
    feature_dat <- FetchData(seurat_obj, feature) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")
    embed_dat <- left_join(embed_dat, feature_dat, by = "cell")
  }

  if (!is.null(plot_dat)){
    embed_dat <- left_join(embed_dat, plot_dat, by = "cell")
  }

  color_aes_str <- feature

  color_aes_str_q <- quo(color_aes_str)
  embed_dat <- embed_dat %>% arrange_at(.vars = color_aes_str)

  p <- ggplot(embed_dat,
              aes_string(xcol, ycol)) +
    geom_point(aes_string(color = paste0("`", color_aes_str, "`")),
               size = pt_size,
               alpha = pt_alpha)

  ## discrete or continuous data?
  if (typeof(embed_dat[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(embed_dat[[feature]])) {
    discrete <- T
  } else {
    discrete <- F
  }

  ## increase legend size
  if (discrete) {
    p <- p + guides(colour = guide_legend(override.aes = list(size = 4))) +
      theme(legend.title = element_blank())
  }

  if (label_text) {
    if(discrete) {
      embed_mean_dat <- embed_dat %>%
        group_by_at(vars(one_of(feature))) %>%
        summarize(med_dim_1 = median(get(xcol)),
                  med_dim_2 = median(get(ycol)))

      p <- p +
        geom_text(data = embed_mean_dat,
                  aes_string(x = "med_dim_1",
                             y = "med_dim_2",
                             label = feature),
                  size = label_size,
                  color = label_color)
    } else {
      warning("label_text not compatible with continuous features")
    }
  }

  ## handle legend limit
  if (is.null(max_y) & !discrete) {
    max_y <- c(0, max(embed_dat[[color_aes_str]]))
  } else if (discrete & is.null(max_y)){
    max_y <- c(NA, NA)
  }

  # loupe-like colors
  cols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])

  #handle legend name
  if(is.null(legend_title)) legend_title <- color_aes_str

  ## handle zero expression
  if (!all(is.na(max_y)) && all(max_y == c(0, 0))){
    p <- p + scale_color_gradient(low = cols[1], high = cols[1], name = legend_title)
    return(p)
  }

  ## handle colors
  if (is.null(.cols) && !discrete){
    if (palette_type == "viridis") {
      p <- p + scale_color_viridis(discrete = F,
                                   direction = -1,
                                   option = col_pal,
                                   limits = max_y, name = legend_title)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1, name = legend_title)
    } else if (palette_type == "cloupe") {
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = cols, name = legend_title)
    }
  } else if (!is.null(.cols) && !discrete){
    p <- p + scale_color_gradientn(limits = max_y,
                                   colors = .cols, name = legend_title)
  } else {

    if(!is.null(.cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = .cols,
        name = legend_title
      )
    } else {
      p <- p + scale_color_manual(
        values = discrete_palette_default,
        name = legend_title
      )
    }
  }
  p
}


plot_umap <- function(...){
  plot_feature(..., embedding = "umap")
}

plot_tsne <- function(...){
  plot_feature(..., embedding = "tsne")
}

plot_pca <- function(...){
  plot_feature(..., embedding = "pca")
}

plot_violin <- function(df, .x, .y,
                        .fill = NULL,
                        .size = 0.50,
                        .width = 1,
                        .scale = "width",
                        .alpha = 1,
                        cols = ggplot2::scale_fill_viridis_d(),
                        single_col = NULL,
                        jitter = F,
                        rotate_x_text = TRUE,
                        arrange_by_fill = TRUE){

  if (arrange_by_fill && !is.null(.fill)){
    tmp <- sym(.fill)
    df <- arrange(df, !!tmp)
    df[[.x]] <- factor(df[[.x]], levels = unique(df[[.x]]))
  }

  p <- ggplot(df, aes_string(x = .x, y = .y))

  if (jitter){
    p <- p  + geom_jitter(size = 0.1, alpha = 0.2, color = "black")
  }

  if (!is.null(single_col)){
    p <- p +
      geom_violin(size = .size,
                  scale = .scale,
                  fill = single_col,
                  alpha = .alpha)
  } else {
    p <- p +
      geom_violin(aes_string(fill = .fill),
                  size = .size,
                  scale = .scale,
                  alpha = .alpha) +
      cols
  }

  if(rotate_x_text){
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.5))
  }
  p <- p + theme(legend.title = element_blank())
  p
}



ilovehue_pal <- c(
   "#e03d6e",
   "#e27c8b",
   "#a64753",
   "#da453e",
   "#db8364",
   "#a54423",
   "#dc652e",
   "#de8c31",
   "#d2a46c",
   "#8f672b",
   "#cea339",
   "#b2b939",
   "#717822",
   "#627037",
   "#a3b46c",
   "#7ba338",
   "#67c042",
   "#3d8829",
   "#35773e",
   "#55c267",
   "#5ca76a",
   "#277257",
   "#5fcea4",
   "#399d82",
   "#40c2d1",
   "#5099cf",
   "#7490df",
   "#615ea5",
   "#716bdf",
   "#c291d6",
   "#984db6",
   "#d558c2",
   "#e17fc0",
   "#995580",
   "#bd3c80"
)
get_distinct_cols <- function(vec, seed = 42) {

  seq_col_pals <- c("Blues", "Greens", "Oranges", "Purples", "Reds", "Greys")
  #seq_cols <- map(seq_col_pals, ~brewer.pal(9, .x) %>% .[1:9] %>% rev(.))

  vec <- sort(vec)
  n_needed <- rle(as.character(vec))$lengths
  n_groups <- length(levels(factor(vec)))

  if(n_groups > 6){
    stop("not enough palettes for ", n_groups, " groups", call. = FALSE)
  }

  seq_col_pals <- seq_col_pals[order(n_needed, decreasing = T)]

  vals <- list()
  for (i in 1:n_groups){
    n <- n_needed[i]
    cols <- suppressWarnings(brewer.pal(n, seq_col_pals[i]))
    if (n < 3){
      cols <- cols[n:3]
    }
    vals[[i]] <- cols
  }
  unlist(vals)
}



set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}


#' Extract out reduced dimensions and cell metadata to tibble
#'
#' @param obj Seurat Object
#' @param embedding dr slot to extract (defaults to all embeddings (2D))
#'
get_metadata <- function(obj, embedding = NULL) {

  mdata <- as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (!embedding %in% names(obj@reductions)) {
      stop(paste0(embedding, " not found in seurat object"), call. = FALSE)
    }

    embed_dat <- obj@reductions[[embedding]]@cell.embeddings %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

  } else {
    embed_dat <- map(names(obj@reductions),
                         ~obj@reductions[[.x]]@cell.embeddings[, 1:2]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

  }

  embed_dat <- left_join(mdata,
                         embed_dat,
                         by = "cell")
  embed_dat
}


#' plot feature across multiple panels split by group
plot_features_split <- function(sobj, feature, group = "orig.ident",
                                embedding = "umap", cols = NULL, add_title = FALSE,
                                ...) {

  # get max coordinates
  dim_reduc <- sobj@reductions[[embedding]]@cell.embeddings[, 1:2]
  x_lims <- c(min(dim_reduc[, 1]), max(dim_reduc[, 1]))
  y_lims <- c(min(dim_reduc[, 2]), max(dim_reduc[, 2]))

  groups <- sort(unique(sobj@meta.data[[group]]))

  if(!is.null(cols)) {
     cols <- cols[1:length(groups)]
     plts <- map2(groups, cols, function(x, y) {
                cells <- rownames(sobj@meta.data)[sobj@meta.data[[group]] == x]
                plot_feature(sobj,
                             feature = feature,
                             embedding = embedding,
                             cell_filter = cells,
                             .cols = y,
                             ...) +
                  coord_cartesian(xlim = x_lims, y = y_lims)
              })
  } else {
    plts <- map(groups, function(x) {
      cells <- rownames(sobj@meta.data)[sobj@meta.data[[group]] == x]
      plot_feature(sobj,
                   feature = feature,
                   embedding = embedding,
                   cell_filter = cells,
                   ...) +
        coord_cartesian(xlim = x_lims, y = y_lims)
    })
  }

  if(add_title){
    plts <- map2(plts, groups, ~.x + labs(title = .y))
  }
  plts
}
