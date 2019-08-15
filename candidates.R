
source("utility.R", local=TRUE)
source("graph_visualization.R")


analyze_components = function(dataset, labels, output_prefix, auroc_threshold = c(0.6)) {    
    is_nonzero_cell = Matrix::colSums(assay(dataset)) > 0
    dataset = dataset[, is_nonzero_cell]
    labels = labels[is_nonzero_cell]
    
    best_hits = compute_best_hits(dataset, labels)
    
    for (threshold in auroc_threshold) {
        output_dir = paste0(output_prefix, threshold*100)
        dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)

        write.table(best_hits, file.path(output_dir, "best_hits.txt"))

        components = extract_components(best_hits, threshold)
        plot_components(best_hits, components$modules, output_dir)
        write_component_summary(components, output_dir)
        export_filtered_best_hits(best_hits, components$modules, output_dir)

        graph = make_directed_graph(best_hits, threshold, 1)
        graph = color_graph(graph, paste(dataset$study_id, labels, sep="|"))
        pdf(file.path(output_dir, "graph_visualization.pdf"))
        plot_directed_graph(graph, 1)
        dev.off()
    }
}

compute_best_hits = function(dataset, labels, one_vs_one = TRUE) {
  normalized_data = normalize_cols(assay(dataset))
  colnames(normalized_data) = paste(dataset$study_id, labels, sep = "|")
  voter_id = design_matrix(colnames(normalized_data))
  voters = normalized_data %*% voter_id
  result = c()
  for (study in unique(dataset$study_id)) {
    candidates = normalized_data[, dataset$study_id == study]
    votes = crossprod(candidates, voters)
    aurocs = compute_aurocs(votes, design_matrix(rownames(votes)))
    if (one_vs_one) {
        result = rbind(result, compute_1v1_aurocs(votes, aurocs))
    } else {
        result = rbind(result, aurocs)
    }
  }
  result = result[, rownames(result)]
  return(result)
}

compute_1v1_aurocs = function(votes, aurocs) {
  result = 0*aurocs
  for (i in seq_len(ncol(aurocs))) {
    top_candidates = find_top_candidates(votes[,i], aurocs[,i])
    result[top_candidates$best, i] = top_candidates$score
    result[top_candidates$second, i] = 1-top_candidates$score
  }
  return(result)
}

find_top_candidates = function(votes, aurocs) {
  candidates = extract_top_candidates(aurocs, 5)
  best = candidates[1]
  votes_best = votes[names(votes) == best]
  score = 1
  second_best = candidates[2]
  for (i in seq(2, length(candidates))) {
    contender = candidates[i]
    votes_contender = votes[names(votes) == contender]
    auroc = c(compute_aurocs(
      as.matrix(c(votes_best, votes_contender)),
      as.matrix(rep(c(1,0), c(length(votes_best), length(votes_contender))))
    ))
    if (auroc > 0.5) {
      if (auroc < score) {
        score = auroc
        second_best = contender
      }
    } else {
      second_best = best
      best = contender
      score = 1 - auroc
      votes_best = votes_contender
    }
  }
  return(list(score = score, best=best, second=second_best))
}

extract_top_candidates = function(aurocs, n = 10) {
  return(names(head(sort(aurocs, decreasing=TRUE), n = n)))
}

extract_components = function(best_hits, threshold = 0) {
  comp = igraph::components(make_graph(best_hits, threshold))
  modules = list()
  outliers = c()
  for (i in seq_len(comp$no)) {
    members = names(which(comp$membership == i))
    if (length(members) > 1) {
      modules[[length(modules)+1]] = members
    } else {
      outliers = c(outliers, members)
    }
  }
  return(list(modules = modules, outliers = outliers))
}

make_graph = function(best_hits, threshold = 0) {
  adj = 0*best_hits
  adj[best_hits > threshold] = 1
  adj = adj * t(adj)
  igraph::graph_from_adjacency_matrix(adj)
}

write_component_summary = function(components, directory = ".") {
  f = file(file.path(directory, "component_summary.txt"), "w")
  writeLines("Label\tComponent", f)
  for (i in seq_along(components$modules)) {
    writeLines(paste(components$modules[[i]], i, sep="\t"), f)
  }
  if (length(components$outliers) > 0) {
    writeLines(paste(components$outliers, -1, sep="\t"), f)
  }
  close(f)
}

plot_components = function(best_hits, components, directory = ".", reorder=FALSE) {
  if (length(components) == 0) { return(); }
             
  dendrogram = if (reorder) "both" else "none"
  color_scale = c(RColorBrewer::brewer.pal(8, "Set2"),
                   RColorBrewer::brewer.pal(9, "Set1"),
                   RColorBrewer::brewer.pal(12, "Set3"))
  breaks = c(0, 0.5, 0.7, 0.9, 0.95, 0.99, 1)
  auroc_cols = colorRampPalette(c("white", "blue"))(length(breaks)-1)
  study_ids = unique(get_study_id(unlist(components)))
  study_cols = color_scale[seq_along(study_ids)]
  names(study_cols) = study_ids
  for (i in seq_along(components)) {
    c = components[[i]]
    dat = best_hits[c, c]
    comp_cols = study_cols[get_study_id(rownames(dat))]
    comp_cell_types = get_cell_type(rownames(dat))
    if (reorder) {
      reorder = as.dendrogram(hclust(as.dist(1-dat), method="average"))
    }
    pdf(file.path(directory, paste0("component_", i, ".pdf")))
    gplots::heatmap.2(
      dat, margins = c(10,10),
      labRow = comp_cell_types, labCol = comp_cell_types,
      key.xlab="AUROC", key.title=NA, cexRow = .7, cexCol = .7,
      trace = "none", col = auroc_cols, breaks = breaks,
      Rowv = reorder, Colv = reorder, dendrogram = dendrogram,
      RowSideColors = rev(comp_cols), ColSideColors = comp_cols,
      revC = TRUE
    )
    par(lend = 1)
    legend("topright", inset = c(0, 0),
           legend = names(study_cols),
           col = study_cols, pt.cex = 1, cex = 0.5, lwd = 10, bty="n")
    dev.off()
  }
}

score_components = function(best_hits, components, full_name = FALSE) {
  name = c()
  n_studies = c()
  score = c()
  for (i in seq_along(components)) {
    c = components[[i]]
    if (full_name) {
        name = c(name, c[1])
    } else {
        name = c(name, get_cell_type(c[1]))
    }
    n_studies = c(n_studies, length(unique(get_study_id(c))))
    score = c(score, mean(best_hits[c, c]))
  }
  return(data.frame(n_studies = n_studies, score = score, row.names = name))
}

export_filtered_best_hits = function(best_hits, components, directory=".") {
  labels = unlist(components)
  filtered_best_hits = matrix(0, nrow=length(labels), ncol=length(labels),
                               dimnames=list(labels, labels))
  for (m in components) {
    filtered_best_hits[m,m] = best_hits[m,m]
  }
  write.table(filtered_best_hits, file.path(directory, 'filtered_best_hits.txt'))
}
