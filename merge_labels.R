
source('candidates2.R')

merge_labels <- function(dataset, labels, threshold = 0) {
  labels <- as.character(labels)
  best_hits <- compute_best_hits(dataset, labels)
  orphans <- find_orphans(best_hits, threshold)
  merge_record <- data.frame()
  while (length(orphans) > 0) {
    print(orphans)
    orphan <- orphans[1]
    merge_candidate <- find_merge_candidate(best_hits, orphan)
    new_labels <- replace_labels(dataset$study_id, labels, c(merge_candidate, orphan))
    new_best_hits <- compute_best_hits(dataset, new_labels)
    new_merge_record <- create_merge_record(labels, best_hits, new_best_hits, orphan, merge_candidate)
    if (new_merge_record$post_merge_consistency / new_merge_record$pre_merge_consistency > 0.9) {
      new_merge_record$merge_status = "success"
      best_hits <- new_best_hits
      labels <- new_labels
    } else {
      new_merge_record$merge_status = "rejected"
      dataset <- dataset[, labels != get_cell_type(orphan)]
      labels <- labels[labels != get_cell_type(orphan)]
      best_hits <- compute_best_hits(dataset, labels)
    }
    merge_record <- rbind(merge_record, new_merge_record)
    orphans <- find_orphans(best_hits, threshold)
  }
  write.table(merge_record, "merge_record.csv")
  return(best_hits)
}

find_orphans <- function(best_hits, threshold = 0) {
  best_hits <- best_hits[best_hits > threshold]
  consistency <- diag(best_hits %*% best_hits)
  incoming_votes <- rowSums(best_hits)
  return(names(sort(consistency[(consistency <= 1) & (incoming_votes <= 1)])))
}

find_merge_candidate <- function(best_hits, orphan) {
  common_hits <- crossprod(best_hits)
  orphan_hits <- common_hits[orphan, get_study_id(colnames(common_hits)) == get_study_id(orphan)]
  orphan_hits <- orphan_hits[names(orphan_hits) != orphan]
  friend <- names(which.max(orphan_hits))
  return(friend)
}

create_merge_record <- function(all_labels, best_hits, new_best_hits, orphan, merge_candidate) {
  consistency <- diag(best_hits %*% best_hits)
  new_consistency <- diag(new_best_hits %*% new_best_hits)
  incoming_votes <- rowSums(best_hits)
  common_votes <- crossprod(best_hits)
  result <- data.frame(
    orphan = orphan, consistency = consistency[orphan], size = sum(all_labels == orphan),
    self_consistency = best_hits[orphan, orphan], incoming_votes = incoming_votes[orphan],
    merge_candidate = merge_candidate, common_votes = common_votes[orphan, merge_candidate],
    pre_merge_consistency = consistency[merge_candidate],
    post_merge_consistency = new_consistency[fuse_labels(c(merge_candidate, orphan))],
    merge_status = NA,
    stringsAsFactors = FALSE
  )
  return(result)
}

replace_labels <- function(study_id, original_labels, merge_candidates) {
  result <- as.character(original_labels)
  study <- get_study_id(merge_candidates[1])
  cell_types <- get_cell_type(merge_candidates)
  new_label <- paste(cell_types, collapse = "+")
  result[(study_id == study) & (result %in% cell_types)] <- new_label
  return(result)
}

fuse_labels <- function(labels_to_merge) {
  study <- get_study_id(labels_to_merge[1])
  cell_types <- get_cell_type(labels_to_merge)
  return(paste(study, paste(cell_types, collapse = "+"), sep = "|"))
}
