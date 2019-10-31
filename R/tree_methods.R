
# Tree methods

is_pruned <- function(xtree, node_n, active = NULL){
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  return(!(node_n %in% active_nN) || with(xtree[[node_n]], !any(children_nN %in% active_nN)))
}

prune_node_old <- function(xtree, node_n, active = NULL) {
  # Pruning a node deactivates recursively any of its children
  # and effectively makes the node itself a leaf
  # Returns list with two members:
  # "active": logical vector marking nodes that are "active" after pruning
  # "pruned": logical vector marking nodes that where "pruned" in the process
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

node_substitution_error <- function(xtree, active = NULL) {
  #NB: For testing purposes only 
  # For every active node returns an absolute value of 
  # a change in error due to a node split

  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)

  return(sapply(seq_along(xtree), function(i, tree, active){
    rv <- with(tree[[i]], {if (i %in% leaves) 0 else sum(best_split$error)})
  }, tree = xtree, active = active))
}

subtree_substitution_error <- function(xtree, active = NULL) {
  # NB: For testing purposes only
  # For every active node "i" returns an absolute value of a change 
  # in error due to substitution of a node "i" for the subtree 
  # with a root at node "i"
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)
  
  info_gain <- sapply(active_nN, function(i, tree, active){
    rv <- with(tree[[i]], {if (i %in% leaves) 0 else sum(best_split$error)})
  }, tree = xtree, active = active)

  sub_error <- vector(mode='numeric', length = length(xtree))
  sub_error[!active] <- 0
  for (i in rev(active_nN)) {
    sub_error[i] <- with(xtree[[i]], {if (i %in% leaves) 0 else info_gain[i] - (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
  }
  sub_error
}

complexity <- function(xtree, active = NULL){
  # For testing purposes only
  # Returns vector of complexity values.
  # "Complexity" value of a node "i" is a number of leaves of a 
  # subtree of 'xtree' with the root at node "i"
  
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)
  
  node_complexity <- rep(0, length = length(xtree))
  for (i in rev(active_nN)) {
    node_complexity[i] <- with(xtree[[i]], {if (i %in% leaves) 1 else (node_complexity[children_nN[1]] + node_complexity[children_nN[2]])})
  }
  node_complexity
}

tree_info <- function(xtree, active = NULL) {

  # "tree_info" defines any nested subtree of the "xtree" 
  # and provides some "helper" information about this subtree
  # that allowes other functions to treet the subtree as "xtree"
  
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)
  
  info_gain <- sapply(active_nN, function(i, tree, active){
    rv <- with(tree[[i]], {if (i %in% leaves) 0 else sum(best_split$error)})
  }, tree = xtree, active = active)
  
  sub_error <- vector(mode='numeric', length = length(xtree))
  sub_error[!active] <- 0
  for (i in rev(active_nN)) {
    sub_error[i] <- with(xtree[[i]], {if (i %in% leaves) 0 else info_gain[i] - (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
  }
  
  node_complexity <- rep(0, length = length(xtree))
  for (i in rev(active_nN)) {
    node_complexity[i] <- with(xtree[[i]], {if (i %in% leaves) 1 else (node_complexity[children_nN[1]] + node_complexity[children_nN[2]])})
  }
  return(list(active = active, node_complexity = node_complexity, substitution_error = sub_error, leaves = leaves))
}

prune_node <- function(xtree, node_n, info = NULL) {
  if (is.null(info)) info <- tree_info(xtree = xtree)
  active <- info$active
  active_nN <- seq_along(xtree)[active]
  if (is_pruned(xtree = xtree, node_n = node_n, active = active)) {
    warn(paste("node", node_n, "has been pruned already"))
    return(list(info = info))
  }
  nodes_to_walk <- node_n
  #one surely can use recursion here
  nactive <- active
  while (length(nodes_to_walk > 0)) {
    curr_node <- nodes_to_walk[1]
    nodes_to_walk <- nodes_to_walk[-1]
    curr_node_children <- xtree[[curr_node]]$children_nN
    delta <- curr_node_children[curr_node_children %in% active_nN]
    if (length(delta) > 0) {
      nactive[delta] <- FALSE
      nodes_to_walk <- c(nodes_to_walk, delta)
    }
  }
  pruned  <- active & (!nactive)
  return(list(active = nactive, pruned = pruned))
  
}