
source("utility.R", local = TRUE)

assess_module_labels <- function(dataset, cell_types, modules, hvg) {
  go_terms <- go_slim(rownames(dataset))
  go_terms$hvg <- hvg
  dataset <- dataset[rownames(dataset) %in% unique(unlist(go_terms)),]
  module_labels <- create_module_labels(dataset$study_id, cell_types, modules)
  return(call_my_metaneighbor(dataset, go_terms, module_labels))
}

call_metaneighbor <- function(dataset, go_terms, module_labels) {
  aurocs <- MetaNeighbor::MetaNeighbor(
    dat = dataset,
    experiment_labels = as.numeric(factor(dataset$study_id)),
    celltype_labels = module_labels,
    genesets = go_terms,
    bplot = FALSE,
    fast_version = TRUE
  )
  return(aurocs)
}

call_my_metaneighbor <- function(dataset, go_terms, module_labels) {
  dat <- assay(dataset)
  result <- list()
  for(l in seq_along(go_terms)){
    print(names(go_terms)[l])
    subdat <- dat[rownames(dat) %in% go_terms[[l]],]
    result[[l]] <- average_auroc(subdat, dataset$study_id, module_labels)
  }
  result <- do.call(rbind, result)
  dimnames(result) <- list(names(go_terms), colnames(module_labels))
  return(result)
}

average_auroc <- function(dat_sub, study_id, celltype_labels) {
  dat_sub <- as.matrix(dat_sub)
  nonzero_cells <- colSums(dat_sub) > 0
  dat_sub <- dat_sub[, nonzero_cells]
  study_id <- study_id[nonzero_cells]
  celltype_labels <- celltype_labels[nonzero_cells,, drop = FALSE]
  normalized_data <- normalize_cols(dat_sub)
  aurocs <- c()
  for (study in unique(study_id)) {
    is_study <- study_id == study
    votes <- crossprod(
      normalized_data[, is_study],
      normalized_data[, !is_study] %*% celltype_labels[!is_study,]
    )
    all_aurocs <- compute_aurocs(votes, celltype_labels[is_study,])
    aurocs <- cbind(aurocs, diag(all_aurocs))
  }
  return(round(rowMeans(aurocs, na.rm = TRUE),3))
}

go_slim <- function(known_genes, min = 10, max = 1000) {
  go_slim <- readRDS("go_slim_ensembl.rds")
  go_slim <- lapply(go_slim, function(g) g[g %in% known_genes])
  group_length <- sapply(go_slim, length)
  return(go_slim[group_length > min & group_length < max])
}

create_module_labels <- function(study_id, cell_types, modules, design=TRUE) {
  result <- paste(study_id, cell_types, sep = "|")
  for (i in seq_along(modules)) {
    m <- modules[[i]]
    result[result %in% m] <- paste("merged", get_cell_type(m[1]), sep = "|")
  }
  if (design) {
    result <- design_matrix(result)
    result <- result[, startsWith(colnames(result), "merged")]
  }
  return(result)
}

go_names <- function() {
  load(file.path('~/data/annotations/GO.human.Rdata'))
  return(voc)
}

go_term_to_name <- function(go_terms, go_name) {
  result <- go_name$V2[match(go_terms, go_name$V1)]
  return(as.character(result))
}

sort_terms <- function(aurocs) {
  result <- lapply(seq_len(ncol(aurocs)), function(i) sort(aurocs[,i], decreasing = TRUE))
  names(result) <- colnames(aurocs)
  return(result)
}

write_top_terms <- function(aurocs, filename) {
  top_terms <- sort_terms(aurocs)
  cat('', file = filename)
  for (i in seq_along(top_terms)) {
    top <- top_terms[[i]][1:5]
    cat('\n', file = filename, append = TRUE)
    cat(names(top_terms)[i], file = filename, sep = '', append = TRUE)
    cat('\n', file = filename, append = TRUE)
    cat(names(top), file = filename, sep = ',', append = TRUE)
    cat('\n', file = filename, append = TRUE)
    cat(top, file = filename, sep = ',', append = TRUE)
    cat('\n', file = filename, append = TRUE)
  }
}
