
library(MetaNeighbor)

source("visualization.R")

iterative_metaneighbor <- function(dataset, depth = 3, output_dir, suffix = "") {
  aurocs <- MetaNeighborUS(
    dat = dataset, var_genes = rownames(dataset),
    study_id = dataset$study_id, cell_type = as.character(dataset$cluster_label),
    fast_version = TRUE
  )
  if ((length(aurocs) < 2) || (any(!is.finite(aurocs)))) { return() }
  pdf(file.path(output_dir, paste0("aurocs", suffix, ".pdf")))
  plot_NV_heatmap(aurocs)
  dev.off()
  if (depth > 1) {
    new_sets <- divide_dataset(dataset, aurocs)
    iterative_metaneighbor(new_sets[[1]], depth-1, output_dir, paste0(suffix, 1))
    iterative_metaneighbor(new_sets[[2]], depth-1, output_dir, paste0(suffix, 2))
  }
}

divide_dataset <- function(dataset, aurocs, cell_labels) {
  data_labels <- paste(dataset$study_id, cell_labels, sep = "|")
  group_labels <- divide_labels(aurocs)
  return(list(dataset[, data_labels %in% group_labels[[1]]],
              dataset[, data_labels %in% group_labels[[2]]]))
}

divide_labels <- function(auroc_table) {
  groups <- cutree(hclust(auroc_as_distance(auroc_table)), 2)
  return(list(names(groups)[groups == 1],
              names(groups)[groups == 2]))
}

auroc_as_distance <- function(auroc_table) {
  # symmetric_auroc <- pmin(auroc_table, t(auroc_table))
  return(as.dist(1-auroc_table))
}
