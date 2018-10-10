
source('utility.R', local=TRUE)

compute_best_hits <- function(dataset, labels) {
  normalized_data <- normalize_cols(assay(dataset))
  colnames(normalized_data) <- paste(dataset$study_id, labels, sep = "|")
  voter_id <- design_matrix(colnames(normalized_data))
  voters <- normalized_data %*% voter_id
  result <- c()
  for (study in unique(dataset$study_id)) {
    candidates <- normalized_data[, dataset$study_id == study]
    votes <- crossprod(candidates, voters)
    aurocs <- compute_aurocs(votes, design_matrix(rownames(votes)))
    #result <- rbind(result, aurocs)
    result <- rbind(result, compute_1v1_aurocs(votes, aurocs))
  }
  result <- result[, rownames(result)]
  return(result)
}

compute_1v1_aurocs <- function(votes, aurocs) {
  result <- 0*aurocs
  for (i in seq_len(ncol(aurocs))) {
    top_candidate <- find_top_candidate(votes[,i], aurocs[,i])
    result[names(top_candidate),i] <- top_candidate
  }
  return(result)
}

find_top_candidate <- function(votes, aurocs) {
  candidates <- extract_top_candidates(aurocs, 5)
  best <- candidates[1]
  votes_best <- votes[names(votes) == best]
  score <- 1
  for (i in seq(2, length(candidates))) {
    contender <- candidates[i]
    votes_contender <- votes[names(votes) == contender]
    auroc <- c(compute_aurocs(
      as.matrix(c(votes_best, votes_contender)),
      as.matrix(rep(c(1,0), c(length(votes_best), length(votes_contender))))
    ))
    if (auroc > 0.5) {
      if (auroc < score) { score <- auroc }
    } else {
      best <- contender
      score <- 1 - auroc
      votes_best <- votes_contender
    }
  }
  names(score) <- best
  return(score)
}

extract_top_candidates <- function(aurocs, n = 10) {
  return(names(head(sort(aurocs, decreasing=TRUE), n = n)))
}

extract_components <- function(best_hits, threshold = 0) {
  comp <- igraph::components(make_graph(best_hits, threshold))
  modules <- list()
  outliers <- c()
  for (i in seq_len(comp$no)) {
    members <- names(which(comp$membership == i))
    if (length(members) > 1) {
      modules[[length(modules)+1]] <- members
    } else {
      outliers <- c(outliers, members)
    }
  }
  return(list(modules = modules, outliers = outliers))
}

make_graph <- function(best_hits, threshold = 0) {
  adj <- 0*best_hits
  adj[best_hits > threshold] <- 1
  adj <- adj * t(adj)
  igraph::graph_from_adjacency_matrix(adj)
}

plot_components <- function(best_hits, components, directory = ".") {
  colfunc <- colorRampPalette(c("white", "blue"))
  breaks <- c(0, 0.5, 0.7, 0.9, 1)
  for (i in seq_along(components)) {
    c <- components[[i]]
    dat <- best_hits[c, c]
    reorder_entries <- as.dendrogram(hclust(as.dist(1-dat), method="average"))
    pdf(file.path(directory, paste0("component_", i, ".pdf")))
    gplots::heatmap.2(
      best_hits[c, c], margins = c(10,10),
      key.xlab="AUROC", key.title=NA, cexRow = .7, cexCol = .7,
      trace = "none", col = colfunc(length(breaks)-1), breaks = breaks,
      Rowv = reorder_entries, Colv = reorder_entries
    )
    dev.off()
  }
}
