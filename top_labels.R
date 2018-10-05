
extract_consistent_subgroups <- function(best_hits) {
  consistency <- sort(diag(best_hits %*% best_hits), decreasing = TRUE)
  modules <- list()
  orphans <- c()
  while (length(consistency) > 0) {
    if (consistency[1] > 1) {
      candidate <- names(consistency)[1]
      neighbors <- names(which(best_hits[, candidate] > 0))
      new_module <- unique(c(candidate, neighbors))
      new_module <- new_module[new_module %in% names(consistency)]
      if (length(new_module) > 1) {
        modules[[length(modules)+1]] <- new_module
      } else {
        orphans <- c(orphans, new_module)
      }
      consistency <- consistency[!(names(consistency) %in% new_module)]
    } else {
      orphans <- c(orphans, names(consistency))
      consistency <- c()
    }
  }
  return(list(modules = modules, orphans = orphans))
}
