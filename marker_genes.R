
library(tidyverse)

source('utility.R')

de_replicability <- function(dataset, study_ids, labels) {
  individual_stats <- list()
  for (study in unique(study_ids)) {
    individual_stats[[study]] <- compute_de_stats(dataset[, study_ids == study], labels[study_ids == study])
  }
  return(list(individual_stats = individual_stats,
              meta_stats = combine_de_stats(individual_stats, labels)))
}

filter_genes <- function(meta_stats, rank_threshold = 0.9, fc_threshold = 4, sparsity_threshold = 0.5) {
  gene_names <- rownames(meta_stats$rank)
  rank_threshold <- rank_threshold * length(gene_names)
  result <- list()
  for (cell_type in colnames(meta_stats$rank)) {
    result[[cell_type]] <- gene_names[
      meta_stats$rank[, cell_type] > rank_threshold &
      meta_stats$fold_change[, cell_type] > fc_threshold &
      meta_stats$expressing_cells[, cell_type] > sparsity_threshold
    ]
  }
  result <- lapply(result, function(r) r[!is.na(r)])
  return(result)
}

compute_de_stats <- function(dataset, cell_types) {
  cell_types <- as.character(cell_types)
  counts <- as.matrix(assay(dataset))
  cpm <- convert_to_cpm(counts)
  return(list(
    aurocs = gene_aurocs(counts, cell_types),
    fold_change = fold_change(cpm, cell_types),
    fold_change_log = fold_change(log1p(cpm), cell_types),
    expressing_cells = expressing_cells(counts, cell_types)
  ))
}

convert_to_cpm <- function(M) {
  return(t(t(M) / colSums(M)) * 1000000)
}

gene_aurocs <- function(dataset, cell_types) {
  result <- compute_aurocs(t(dataset), design_matrix(cell_types))
  return(t(result))
}

fold_change <- function(dataset, cell_types) {
  centroids <- as.matrix(dataset) %*% design_matrix(cell_types, scale = TRUE)
  result <- (ncol(centroids)-1) * (centroids / (rowSums(centroids) - centroids))
  return(result)
}

expressing_cells <- function(dataset, cell_types) {
  result <- as.matrix(dataset > 0) %*% design_matrix(cell_types, scale = TRUE)
  return(result)
}

combine_de_stats <- function(list_datasets, cell_types) {
  label_levels <- unique(cell_types)
  genes <- rownames(list_datasets[[1]][[1]])
  aurocs <- matrix(0, nrow = length(genes), ncol = length(label_levels),
                   dimnames = list(genes, label_levels))
  fold_change <- expressing_cells <- aurocs
  n_studies <- rep(0, length(label_levels))
  names(n_studies) <- label_levels
  for (dataset in list_datasets) {
    study_labels <- colnames(dataset$aurocs)
    aurocs[, study_labels] <- aurocs[, study_labels] + apply(dataset$aurocs, 2, rank)
    fold_change[, study_labels] <- fold_change[, study_labels] + dataset$fold_change
    expressing_cells[, study_labels] <- expressing_cells[, study_labels] + dataset$expressing_cells
    n_studies[study_labels] <- n_studies[study_labels] + 1
  }
  return(list(n_studies = n_studies,
              rank = sweep(aurocs, 2, n_studies, "/"),
              fold_change = sweep(fold_change, 2, n_studies, "/"),
              expressing_cells = sweep(expressing_cells, 2, n_studies, "/")))
}

plot_heatmap <- function(gene_expression) {
  stat_max <- max(abs(min(gene_expression)), abs(max(gene_expression)))
  breaks <- seq(-stat_max, stat_max, length = 21)
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20))
  heatmap.2(gene_expression, Rowv = TRUE, Colv = FALSE, dendrogram = 'none', trace = 'none', breaks = breaks, col = cols)
}

plot_violins <- function(sce_data, markers, labels) {
  result <- ggplot(gather_expression(sce_data, markers, labels)) +
            geom_violin(aes(x = gene, y = log1p(expression))) +
            facet_grid(label ~ .) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            theme(strip.text.y = element_text(angle = 0))
  print(result)
}

gather_expression <- function(sce_data, markers, labels) {
  result <- list()
  labels <- as.factor(labels)
  for (l in levels(labels)) {
    result[[length(result) + 1]] <- data.frame(
      gene = factor(markers, levels = markers),
      label = l,
      expression = c(assay(sce_data[, labels == l]))
    )
  }
  return(Reduce(rbind, result))
}

plot_dots <- function(expression, markers, labels) {
  average_expression <- as.data.frame(expression %*% design_matrix(labels, scale = TRUE))
  average_expression$gene <- factor(markers, levels = markers)
  average_expression <- average_expression %>% gather(label, average_expression, -gene)
  expressing <- as.data.frame(expressing_cells(expression, labels))
  expressing$gene <- factor(markers, levels = markers)
  expressing <- expressing %>% gather(label, expressing_cells, -gene)
  to_plot <- average_expression %>% left_join(expressing, by = c('gene', 'label'))
  result <- ggplot(to_plot) +
    geom_point(aes(x = label, y = gene, col = average_expression, size = expressing_cells)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(result)
}

gene_sets_to_predictors <- function(gene_sets, dataset) {
  return(crossprod(sapply(gene_sets, function(g) rownames(dataset) %in% g), assay(dataset)))
}

extract_top_genes <- function(aurocs, n = 5, decreasing = TRUE) {
  result <- list()
  for (i in 1:nrow(aurocs)) {
    result[[i]] <- top_labels(aurocs[i,], n, decreasing = decreasing)
  }
  names(result) <- rownames(aurocs)
  return(result)
}

top_labels <- function(values, n = 5, decreasing = TRUE) {
  return(names(head(sort(values, decreasing = decreasing), n)))
}
