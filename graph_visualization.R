
make_directed_graph = function(best_hits, low_threshold = 0, high_threshold = 1) {
    filtered_hits = best_hits
    filtered_hits[best_hits>high_threshold | best_hits < low_threshold] = 0
    result = igraph::graph_from_adjacency_matrix(t(filtered_hits), weighted = TRUE)
    result = igraph::simplify(result, remove.loops = TRUE)
    return(result)
}

color_graph = function(graph, full_label_list) {
    vertex_colors = make_vertex_colors(graph)
    vertex_study = factor(get_study_id(igraph::V(graph)$name))
    levels(vertex_study) = vertex_colors[levels(vertex_study)]
    
    cluster_size = table(full_label_list)
    cluster_size_class = as.numeric(cut(cluster_size, breaks = c(0, 10, 100, max(cluster_size)+1)))
    names(cluster_size_class) = names(cluster_size)

    igraph::V(graph)$color = as.character(vertex_study)
    igraph::V(graph)$label.color = "black"
    igraph::V(graph)$label = get_cell_type(igraph::V(graph)$name)
    igraph::V(graph)$size = cluster_size_class[igraph::V(graph)$name]

    igraph::E(graph)$width = igraph::E(graph)$weight
    igraph::E(graph)$color = c("orange","darkgray")[as.numeric(igraph::E(graph)$weight >= 0.5) + 1]

    return(graph)
}

make_vertex_colors = function(graph) {
    study_ids = unique(get_study_id(igraph::V(graph)$name))
    result = RColorBrewer::brewer.pal(n = 8, "Set2")[seq_along(study_ids)]
    names(result) = study_ids
    return(result)
}

plot_directed_graph = function(graph, size_factor=1) {
    vertex_colors = make_vertex_colors(graph)
    vertex_size = igraph::V(graph)$size * size_factor
    edge_width = igraph::E(graph)$width * size_factor
    plot(graph, vertex.label.cex=0.2*size_factor, edge.arrow.size=.05*size_factor**1.5, vertex.frame.color = NA,
         vertex.label.font=2, vertex.size = vertex_size, edge.width = edge_width)
    legend("topright", legend = names(vertex_colors), pt.bg = vertex_colors, pt.cex = 1, cex = 0.5, bty="n", pch=21)
}

extend_set_to_neighbors = function(graph, initial_set, max_neighbor_distance=2) {
    A = igraph::as_adj(graph)
    A[A>0] = 1
    diag(A) = 1
    V = as.numeric(igraph::V(graph)$name %in% initial_set)
    result = V
    for (i in seq_len(max_neighbor_distance)) {
        result = crossprod(A, result)
    }
    return(rownames(result)[as.logical(result)])
}

make_subgraph = function(graph, vertices) {
    return(igraph::induced_subgraph(graph, vertices))
}
