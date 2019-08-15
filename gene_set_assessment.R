
source("utility.R", local = TRUE)
current_dir__ = getwd()
common_dir__ = "~/projects/common"
setwd(common_dir__)
source("go_annotation.R")
source("hgnc.R")
source("identifier_conversion.R")
setwd(current_dir__)

assess_module_labels = function(dataset, cell_types, modules, hvg) {
  go_terms = get_filtered_go_mouse(rownames(dataset))
  go_terms$hvg = hvg
  dataset = dataset[rownames(dataset) %in% unique(unlist(go_terms)),]
  module_labels = create_module_labels(dataset$study_id, cell_types, modules)
  return(call_my_metaneighbor(dataset, go_terms, module_labels))
}

get_filtered_go_mouse = function(known_genes, min_size = 20, max_size = 1000) {
    current_dir__ = getwd()
    result = go_mouse()
    result = result[rownames(result) %in% known_genes,]
    set_size = Matrix::colSums(result)
    result = result[, set_size > min_size & set_size < max_size]
    result = apply(result, 2, function(x) names(which(x)))
    return(result)
}

call_metaneighbor = function(dataset, go_terms, module_labels) {
  aurocs = MetaNeighbor::MetaNeighbor(
    dat = dataset,
    experiment_labels = as.numeric(factor(dataset$study_id)),
    celltype_labels = module_labels,
    genesets = go_terms,
    bplot = FALSE,
    fast_version = TRUE
  )
  return(aurocs)
}

call_my_metaneighbor = function(dataset, go_terms, module_labels) {
  dat = assay(dataset)
  result = list()
  for(l in seq_along(go_terms)){
    print(names(go_terms)[l])
    subdat = dat[rownames(dat) %in% go_terms[[l]],]
    print(dim(subdat))
    result[[l]] = average_auroc(subdat, dataset$study_id, module_labels)
  }
  result = do.call(rbind, result)
  dimnames(result) = list(names(go_terms), colnames(module_labels))
  return(result)
}

average_auroc = function(dat_sub, study_id, celltype_labels) {
  dat_sub = as.matrix(dat_sub)
  nonzero_cells = colSums(dat_sub) > 0
  dat_sub = dat_sub[, nonzero_cells]
  study_id = study_id[nonzero_cells]
  celltype_labels = celltype_labels[nonzero_cells,, drop = FALSE]
  normalized_data = normalize_cols(dat_sub)
  aurocs = c()
  for (study in unique(study_id)) {
    is_study = study_id == study
    votes = crossprod(
      normalized_data[, is_study],
      normalized_data[, !is_study] %*% celltype_labels[!is_study,]
    )
    all_aurocs = compute_aurocs(votes, celltype_labels[is_study,])
    aurocs = cbind(aurocs, diag(all_aurocs))
  }
  return(round(rowMeans(aurocs, na.rm = TRUE),3))
}

go_slim = function(known_genes, min = 10, max = 1000) {
  go_slim = readRDS("go_slim_ensembl.rds")
  go_slim = lapply(go_slim, function(g) g[g %in% known_genes])
  group_length = sapply(go_slim, length)
  return(go_slim[group_length > min & group_length < max])
}

create_module_labels = function(study_id, cell_types, modules, design=TRUE) {
  result = paste(study_id, cell_types, sep = "|")
  for (i in seq_along(modules)) {
    m = modules[[i]]
    result[result %in% m] = paste("merged", get_cell_type(m[1]), sep = "|")
  }
  if (design) {
    result = design_matrix(result)
    result = result[, startsWith(colnames(result), "merged")]
  }
  return(result)
}

assess_module_labels_hgnc = function(dataset, cell_types, modules, hvg) {
  sets = get_filtered_hgnc(rownames(dataset))
  sets$hvg = hvg
  dataset = dataset[rownames(dataset) %in% unique(unlist(sets)),]
  module_labels = create_module_labels(dataset$study_id, cell_types, modules)
  return(call_my_metaneighbor(dataset, sets, module_labels))
}

get_filtered_hgnc = function(known_genes, min_size = 5, max_size = 500) {
    result = hgnc_families()
    rownames(result) = convert_hgnc_to_mgi_symbols(rownames(result))
    result = result[rownames(result) %in% known_genes,]
    set_size = Matrix::colSums(result)
    result = result[, set_size > min_size & set_size < max_size]
    result = apply(result, 2, function(x) names(which(x)))
    return(result)
}

assess_module_labels_josh = function(dataset, cell_types, modules, hvg) {
  sets = get_filtered_josh(rownames(dataset))
  sets$hvg = hvg
  dataset = dataset[rownames(dataset) %in% unique(unlist(sets)),]
  module_labels = create_module_labels(dataset$study_id, cell_types, modules)
  return(call_my_metaneighbor(dataset, sets, module_labels))
}

get_filtered_josh = function(known_genes, min_size = 5, max_size = 500) {
    josh_lists = read.table("josh_gene_lists.txt", stringsAsFactors = FALSE)
    result = split(josh_lists$gene, josh_lists$Gene_family)
    set_size = sapply(result, function(i) sum(i %in% known_genes))
    result = result[set_size > min_size & set_size < max_size]
    return(result)
}
