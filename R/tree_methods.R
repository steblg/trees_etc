
# Tree methods

is_pruned <- function(xtree, node_n, active = NULL){
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  return(!(node_n %in% active_nN) || with(xtree[[node_n]], !any(children_nN %in% active_nN)))
}

prune_node <- function(xtree, node_n, active = NULL) {
  # Pruning a node deactivates recursively any children on this node
  # and effectively makes it a leaf
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  if (is_pruned(xtree = xtree, node_n = node_n, active = active)) {
    warn(paste("node", node_n, "has been pruned already"))
    rv <- list(active = active)
  } else {
    nodes_to_check <- node_n
    #one surely can use recursion here
    nactive <- active
    while (length(nodes_to_check > 0)) {
      curr_node <- nodes_to_check[1]
      nodes_to_check <- nodes_to_check[-1]
      curr_node_children <- xtree[[curr_node]]$children_nN
      delta <- curr_node_children[curr_node_children %in% active_nN]
      if (length(delta) > 0) {
        nactive[delta] <- FALSE
        nodes_to_check <- c(nodes_to_check, delta)
      }
    }
    pruned  <- active & (!nactive)
    return(list(active = nactive, pruned = pruned))
  }
}

tree_leaves <- function(xtree, active  = NULL, index = TRUE){
  # If isTRUE(index) returns index, i.e. nodes numbers, otherwise returns selector, i.e. logical vector of length(xtree)
  if (is.null(active)) 
    leaves_F <- seq_along(xtree)[sapply(xtree, function(x) is.null(x$children_nN))]
  else {
    stopifnot(length(xtree) == length(active))
    active_nN <- seq_along(xtree)[active]
    leaves_F <- active_nN[sapply(active_nN, function(i){with(xtree[[i]], is.null(children_nN) || !any(active[children_nN]))})]
  }
  rv <- leaves_F
  if (!isTRUE(index)) {
    rv <- rep(FALSE, length(xtree))
    rv[leaves_F] <- TRUE
  }
  rv
}

information_gain <- function(xtree, active = NULL) {
  # if (is.null(active)){
    # return(sapply(xtree, function(x) {
      # rv <- with(x, {if (is.null(best_split)) 0 else sum(best_split$error)})
    # }))
  # }
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  return(sapply(seq_along(xtree), function(i, tree, active){
    rv <- with(tree[[i]], {if (is.null(best_split) || !active[i]) 0 else sum(best_split$error)})
  }, tree = xtree, active = active))
}

substitution_error <- function(xtree, active = NULL) {

  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  
  info_gain <- information_gain(xtree, active = active)
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)
  
  sub_error <- vector(mode='numeric', length = length(xtree))
  sub_error[!active] <- 0
  for (i in rev(active_nN)) {
    sub_error[i] <- with(xtree[[i]], {if (i %in% leaves) 0 else info_gain[i] - (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
  }
  sub_error
}

tree_info <- function(xtree, active = NULL) {
  # xtree: is a basic regression/classification tree
  # active_nodes: optional logical vector specifying "active" nodes
  # Nodes becaome inactive as a result of pruning

  # "tree_info" for every node returns information about
  # leaves: TRUE if node has "active" children
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  
  info_gain <- information_gain(xtree, active = active)
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)
  
  sub_error <- vector(mode='numeric', length = length(xtree))
  sub_error[!active] <- 0
  for (i in rev(active_nN)) {
    sub_error[i] <- with(xtree[[i]], {if (i %in% leaves) 0 else info_gain[i] - (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
  }
  return(list(active = active, information_gain = info_gain, substitution_error = sub_error, leaves = leaves))
}

prune_tree <- function(xtree, node_n, tree_info = NULL) {
}