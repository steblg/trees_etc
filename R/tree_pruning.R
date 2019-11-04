# Tree methods

is_pruned <- function(xtree, node_n, active = NULL){
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  return(!(node_n %in% active_nN) || with(xtree[[node_n]], !any(children_nN %in% active_nN)))
}

split_tree <- function(xtree, node_n, active = NULL) {
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
    rv <- with(tree[[i]], {if (i %in% leaves) 0 else error - sum(best_split$error)})
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
    rv <- with(tree[[i]], {if (i %in% leaves) 0 else error - sum(best_split$error)})
  }, tree = xtree, active = active)
  
  # info_gain <- node_substitution_error(xtree = xtree, active = active)

  sub_error <- vector(mode='numeric', length = length(xtree))
  sub_error[!active] <- 0
  for (i in rev(active_nN)) {
    sub_error[i] <- with(xtree[[i]], {if (i %in% leaves) 0 else info_gain[i] + (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
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
  
  info_gain <- sub_error <- node_complexity <- rep(0, length(xtree))
  
  for (i in rev(active_nN)) {
    curr_node <- xtree[[i]]
    info_gain[i] <- with(curr_node, {if (i %in% leaves) 0 else error - sum(best_split$error)}) 
    sub_error[i] <- with(curr_node, {if (i %in% leaves) 0 else info_gain[i] + (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
    node_complexity[i] <- with(curr_node, {if (i %in% leaves) 1 else (node_complexity[children_nN[1]] + node_complexity[children_nN[2]])})
  }
  
  tmp <- rep(FALSE, length(xtree))
  tmp[leaves] <- TRUE
  return(list(active = active, split_info_gain = info_gain, node_complexity = node_complexity, substitution_error = sub_error, leaves = tmp))
}

walk_down <- function(xtree, node_n, active = NULL, index = TRUE){
  # Returns vector of all node indices in the subtree
  # of the "xtree" with a root at "node_n"
  # Return vector does not include "node_n"
  
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  
  nodes_to_walk <- node_n
  #one surely can use recursion here
  nactive <- active
  
  # walking down the subtree that is being pruned
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
  if (isTRUE(index)) pruned <- seq_len(length(xtree))[pruned]
  return(pruned)
}

walk_up <- function(xtree, node_n) {
  # Returns vector of indices of nodes on the path
  # from the "node_n" up to the root node of the "xtree"
  
  tree_path <- numeric(0)
  curr_node <- node_n
  while (!is.null(curr_node)){
    tree_path <- c(tree_path, curr_node)
    curr_node <- xtree[[curr_node]]$parent_nN 
  }
  return(tree_path)
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
  nactive <- active
  
  complexity_delta <- with(info, 1 - node_complexity[node_n])
  substitution_error_delta <- with(info, substitution_error[node_n])

  # walking down the subtree that is being pruned
  pruned <- walk_down(xtree = xtree, node_n = node_n, active = active, index = TRUE)
  info <- within(info, {
    active[pruned] <- FALSE
    leaves[pruned] <- FALSE
    leaves[node_n] <- TRUE
    split_info_gain[pruned] <- 0
    split_info_gain[node_n] <- 0
    node_complexity[pruned] <- 0
    substitution_error[pruned] <- 0
  })
  
  # propagate change up the tree from the node that is being pruned 
  walk_up_path <- walk_up(xtree = xtree, node_n = node_n)
  for ( node_i in walk_up_path){
    info$node_complexity[node_i] <- info$node_complexity[node_i] + complexity_delta
    info$substitution_error[node_i] <- info$substitution_error[node_i] - substitution_error_delta
  }
  
  return(info)
}