library(RColorBrewer)
library(ROCR)

plot_roc <- function(predictors, labels) {
  colors = RColorBrewer::brewer.pal(ncol(predictors), 'Set1')
  for (i in seq_len(ncol(predictors))) {
    pred <- ROCR::prediction(predictors[,i], labels)
    plot(ROCR::performance(pred, 'tpr', 'fpr'), col = colors[i], add = i!=1)
  }
  legend('bottomright', legend = colnames(predictors), fill = colors)
}

#' predictors is a matrix where each column is a predictor and each row is a sample.
#' label_matrix is a binary matrix where columns are labels and each row is a sample.
#' 1 indicates that the sample on this row belongs to the label on this column.
compute_aurocs <- function(predictors, label_matrix) {
  n_positives <- colSums(label_matrix)
  n_negatives <- nrow(label_matrix) - n_positives
  sum_of_positive_ranks <- crossprod(
    label_matrix,
    matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
  )
  colnames(sum_of_positive_ranks) <- colnames(predictors)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

design_matrix <- function(cell_type, scale = FALSE) {
  factors <- levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result <- model.matrix(~cell_type-1)
  } else {
    result <- matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) <- factors
  if (scale) {
    result <- t(t(result) / colSums(result))
  }
  return(result)
}

normalize_cols <- function(M, ranked = TRUE) {
  M <- as.matrix(M)
  if (ranked) {
    M <- matrixStats::colRanks(M, ties.method = "average", preserveShape = TRUE)
  }
  return(scale_cols(M))
}

Rcpp::cppFunction('NumericMatrix scale_cols(NumericMatrix M) {
  NumericMatrix result(M.nrow(), M.ncol());
  for (int j = 0; j < M.ncol(); j++) {
    double m = 0;
    for (int i = 0; i < M.nrow(); i++) { m += M(i,j); }
    m /= M.nrow();
    for (int i = 0; i < M.nrow(); i++) { result(i,j) = M(i,j) - m; }
    double s = 0;
    for (int i = 0; i < M.nrow(); i++) { s += result(i,j) * result(i,j); }
    s = 1 / sqrt(s);
    for (int i = 0; i < M.nrow(); i++) { result(i,j) *= s; }
  }
  return result;
}')

get_study_id <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed=TRUE), "[", 1))
}

get_cell_type <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed=TRUE), "[", -1))
}
