
source("utility.R", local = TRUE)

assess_module_labels <- function(dataset, cell_types, modules, hvg) {
  go_terms <- go_slim(rownames(dataset))
  go_terms$hvg <- hvg
  dataset <- dataset[unique(unlist(go_terms)),]
  aurocs <- MetaNeighbor::MetaNeighbor(
    dat = dataset,
    experiment_labels = as.numeric(factor(dataset$study_id)),
    celltype_labels = create_module_labels(dataset$study_id, cell_types, modules),
    genesets = go_terms,
    bplot = FALSE,
    fast_version = TRUE
  )
  return(aurocs)
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
