
library(RColorBrewer)

extract_label_to_subclass_mapping <- function(dataset) {
  mapping <- colData(dataset)[, c("cluster_label", "subclass_label")]
  mapping <- result[match(unique(result$cluster_label), result$cluster_label), ]
  result <- mapping$subclass_label
  names(result) <- mapping$cluster_label
  return(result)
}

replace_results_labels <- function(mn_results, label_mapping) {
  rownames(mn_results) <- replace_labels(rownames(mn_results), label_mapping)
  colnames(mn_results) <- replace_labels(colnames(mn_results), label_mapping)
  return(mn_results)
}

replace_labels <- function(labels, label_mapping) {
  return(paste(get_study_id(labels), label_mapping[get_cell_type(labels)], sep = "|"))
}

plot_NV_heatmap <- function(
  dat, reorder_entries = TRUE, breaks = seq(0, 1, length = 21), norm = ""
) {
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20))
  row_colors = get_colors(get_cell_type(colnames(dat)))
  col_colors = get_colors(get_study_id(colnames(dat)))
  if (reorder_entries) {
    reorder_entries <- as.dendrogram(hclust(as.dist(1-dat), method="average"))
  }
  if (norm == "rank") {
     dat <- rank_normalize(dat)
  } else if (norm == "log") {
     dat <- log_normalize(dat)
  }
  gplots::heatmap.2(
    dat, margins = c(11,11),
    key = TRUE, keysize = 1, key.xlab="AUROC", key.title="NULL",
    labRow = NA, offsetRow=0.1, labCol = NA, offsetCol=0.1,
    trace = "none", density.info = "none", col = cols, breaks = breaks,
    Rowv = reorder_entries, Colv = reorder_entries, dendrogram = "row",
    RowSideColors = row_colors$colors, ColSideColors = col_colors$colors
  )
  par(lend = 1)
  legend("topright", inset = c(0, .2),
         legend = col_colors$legend,
         col = col_colors$color_scale, pt.cex = 1, cex = 0.5, lwd = 10, bty="n")
  legend("topright", inset = c(0, .3),
         legend = row_colors$legend,
         col = row_colors$color_scale, cex = 0.5, lwd = 10, bty="n")
}

get_colors <- function(labels) {
  labels <- as.factor(labels)
  numeric_labels <- as.numeric(labels)
  color_scale <- c(RColorBrewer::brewer.pal(8, "Set2"),
                   RColorBrewer::brewer.pal(9, "Set1"),
                   RColorBrewer::brewer.pal(12, "Set3"))
  return(list(colors = color_scale[numeric_labels],
              legend = levels(labels),
              color_scale = color_scale[seq_along(levels(labels))]))
}

rank_normalize <- function(matrix_) {
  row_labels <- rownames(result)
  col_labels <- colnames(result)
  result <- matrix(rank(matrix_), nrow = nrow(matrix_))
  rownames(result) <- row_labels
  colnames(result) <- col_labels
  return(result / max(result))
}

log_normalize <- function(matrix_) {
  result <- matrix_
  high <- result[result <= 0.5]
  low <- result[result > 0.5]
  result[result <= 0.5] <- log(high * (1 - high))
  result[result > 0.5] <- -log(low * (1 - low))
  result <- result - min(result)
  return(result / max(result))
}


bplot <- function(nv_mat, hvg=NULL, cex=1) {
  Celltype = rep(colnames(nv_mat),each=dim(nv_mat)[1])
  ROCValues = unlist(lapply(seq_len(dim(nv_mat)[2]), function(i) nv_mat[,i]))
  beanplot::beanplot(
    ROCValues ~ Celltype, border="NA", col="gray", ylab="AUROC", log = "",
    what=c(0,1,1,1), frame.plot = FALSE, las = 2, cex.axis = cex
  )
  if (!is.null(hvg)) {
    points(hvg ~ factor(names(hvg)), pch=17, col="red")
  }
}
