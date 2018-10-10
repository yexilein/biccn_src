
library(SingleCellExperiment)

create_allen_paul <- function() {
  paul <- readRDS('~/data/brain/paul.rds')
  paul_cluster <- sapply(strsplit(colnames(paul), "\\."), head, 1)
  paul_subclass <- paul_cluster
  paul_subclass[paul_subclass %in% c('PVC', 'CHC1', 'CHC2')] <- 'Pvalb'
  paul_subclass[paul_subclass %in% c('MNC', 'LPC')] <- 'Sst'
  paul_subclass[paul_subclass %in% c('ISC', 'CCKC')] <- 'Vip'
  paul_col_data <- data.frame(study_id = 'paul', class_label = 'GABAergic',
                              subclass_label = paul_subclass,
                              cluster_label = paul_cluster)
  result <- add_to_allen(paul, paul_col_data)
  return(result)
}

add_to_allen <- function(other, other_col_data) {
  allen <- readRDS('allen.rds')
  allen_col_data <- colData(allen)[,colnames(other_col_data)]
  genes <- intersect(rownames(allen), rowData(other)$ensembl_gene_id)
  sub_allen <- allen[genes,]
  sub_other <- other[match(genes, rowData(other)$ensembl_gene_id),]
  assay(sub_other) <- Matrix::Matrix(assay(sub_other), sparse=TRUE)
  result <- SingleCellExperiment(
    list(counts = cbind(assay(sub_allen), assay(sub_other))),
    colData = rbind(allen_col_data, other_col_data)
  )
  exprs(result) <- log1p(scater::calculateCPM(result, use_size_factors = FALSE))
  return(result)
}

create_tasic <- function() {
  tasic <- readRDS('~/data/brain/tasic.rds')
  col_data <- data.frame(
    study_id = 'tasic', class_label = tasic$broad_type,
    subclass_label = tasic$broad_type, cluster_label = tasic$primary_type
  )
  levels(col_data$class_label)[levels(col_data$class_label) == 'Unclassified'] <- 'Noise'
  return(add_to_allen(tasic, col_data))
}

create_zeisel <- function() {
  zeisel <- readRDS('~/data/brain/zeisel.rds')
  class_label <- zeisel$level1class
  class_label <- plyr::revalue(class_label, c(
    "astrocytes_ependymal"="Non-neuronal",
    "endothelial-mural"="Non-neuronal",
    "interneurons"="GABAergic",
    "microglia"="Non-neuronal",
    "oligodendrocytes"="Non-neuronal",
    "pyramidal CA1"="Glutamatergic",
    "pyramidal SS"="Glutamatergic"
  ))
  col_data <- data.frame(
    study_id = 'zeisel', class_label = class_label,
    subclass_label = zeisel$level1class, cluster_label = zeisel$level2class
  )
  rownames(col_data) <- colnames(zeisel)
  return(add_to_allen(zeisel, col_data))
}


variable_genes <- function(dataset, sample_size = 50000, i = 1) {
  if (sample_size < ncol(dataset)) {
    subset <- dataset[,sample.int(ncol(dataset), sample_size)]
  } else {
    subset = dataset
  }
  return(MetaNeighbor::variableGenes(subset, i, exp_labels=subset$study_id))
}

scran_variable_genes <- function(dataset, n = 200) {
  if (length(unique(dataset$study_id)) > 1) {
    design <- model.matrix(~as.character(dataset$study_id))
    fit <- scran::trendVar(dataset, design=design, use.spikes=FALSE)
  } else {
    fit <- scran::trendVar(dataset, use.spikes=FALSE)
  }
  result <- scran::decomposeVar(dataset, fit)
  result <- result[order(result$bio, decreasing=TRUE), ]
  return(rownames(result)[seq_len(n)])
}
