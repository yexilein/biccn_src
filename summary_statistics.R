
summary_score = function(best_hits) {
    return(rowSums(best_hits * t(best_hits)))
}

filter_best_hits = function(best_hits) {
    best_hits[best_hits < 0.6] = 0
    best_hits[best_hits > 0] = 1
    return(best_hits * t(best_hits))
}

get_dataset_matrix = function(best_hits) {
    datasets = factor(sapply(strsplit(rownames(best_hits), split = "|", fixed=TRUE), "[", 1))
    dataset_matrix = model.matrix(~datasets-1)
    colnames(dataset_matrix) = levels(datasets)
    return(dataset_matrix)
}

get_dataset_mapping = function(best_hits) {
    result = filter_best_hits(best_hits) %*% get_dataset_matrix(best_hits)
    result[result > 0] = 1
    return(result)
}

plot_dataset_mapping = function(dataset_mapping) {
    UpSetR::upset(as.data.frame(dataset_mapping), order.by = "freq", sets = colnames(dataset_mapping))
}
