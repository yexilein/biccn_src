
source('~/data/annotations/annotations.R')
source('~/data/biccn/parsed_data/zeng.R')

DATA_DIR = '~/data/biccn/parsed_data'

read_data <- function() {
  return(list(
    zeng_cells = readRDS(file.path(DATA_DIR, 'zeng_10x_cells.rds'))
  ))
}

ribosomal_genes <- function() {
  go_notes <- go_mouse()
  mart_notes <- mart_mouse()
  result <- get_GO_genes(go_notes, 'GO:0022626')
  result <- symbol_to_ensembl(result, mart_notes, human=FALSE)$ensembl_gene_id
  return(result)
}

plot_ribosome_expression <- function(full_dataset, ribo_dataset) {
  expressed_ribosomal_proteins <- colSums(ribo_dataset > 0)
  expressed_genes <- colSums(full_dataset > 0)
  plot(expressed_ribosomal_proteins, expressed_genes)
}
