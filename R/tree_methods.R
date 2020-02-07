
      # errorVal <- c(errorVal, sum(sapply(rtree_, function(x){
        # if(x$status == 'P') return(0) else return(x$error)
      # })))


#To be reworked when using S3

node_defn <- function(node, as_str = FALSE){
  slL <- list()
  for(e in node$split_conditions_list){
    slL[[e$feature]] <- c(slL[[e$feature]], list(e))
  }
  collapse_feature <- function(cL){ # elements of cL have the same 'feature'
      format_rule <- function(x) {
        if (!as_str) return(x)
        if (length(x) == 1) {
          rv <- with(x[[1]], paste("(", feature, operation, split_value, ")"))
        } else if (length(x) == 2 && x[[1]]$operation %in% c("<", ">=")) {
          if (x[[1]]$operation == "<")
            rv <- paste("(", x[[2]]$split, "<=", x[[1]]$feature, x[[1]]$operation, x[[1]]$split, ")")
          else
            rv <- paste("(", x[[1]]$split, "<=", x[[1]]$feature, x[[2]]$operation, x[[2]]$split, ")")  
        }
        rv
      }
      if(length(cL) == 1) return(format_rule(cL))
      opF <- factor(sapply(cL, function(x) x$operation))
      L <- split(cL, opF)
      names(L) <- NULL
      collapse_F <- function(y) { # elements of y have same 'feature' and the same 'operation'
        if(y[[1]]$operation %in% c("==", '!=')) return(y)
        if(y[[1]]$operation == "<") return(list(feature = y[[1]]$feature, operation = y[[1]]$operation, split_value = min(sapply(y, function(x) x$split_value))))
        if(y[[1]]$operation == ">=") return(list(feature = y[[1]]$feature, operation = y[[1]]$operation, split_value = max(sapply(y, function(x) x$split_value))))
        stop(paste(paste("'", y[[1]]$operation, "'", sep =""), "is not supported"))
      }
      rv <- format_rule(lapply(L, collapse_F))
      return(rv)
  }
  rv <- lapply(slL, collapse_feature)
  if (isTRUE(as_str)) rv <- paste( unlist(rv), collapse = " & ")
  rv
}

# node_defn_str <- function(node, ...){
  # definition <- unlist(node_defn(node = node, ...), recursive = FALSE, use.names = FALSE)
  # return(paste(sapply(definition, function(x) with(x, paste("(", feature, operation, split_value, ")"))), collapse = " & "))
# }
    
to_data_tree_node <- function(rtree_node, abbrev = TRUE) {
  library(data.tree)

  if(isTRUE(abbrev)){
    nm <- with(rtree_node, split_conditions[length(split_conditions)])
  } else {
    nm <- node_defn(rtree_node, as_str = TRUE)
  }
  rv <- Node$new(nm)
  for (nm in names(rtree_node)) {
    rv[[nm]] <- rtree_node[[nm]]
  }
  rv
}

to_data_tree <- function(xtree, active = NULL){
  stopifnot(is.null(active) || length(active) == length(xtree))
  library(data.tree)
  
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  flat_tree <- lapply(xtree, to_data_tree_node)
  assembled_tree <- flat_tree[[1]]
  for(i in 2 : length(flat_tree)){
    if (active[i]) flat_tree[[flat_tree[[i]]$parent_nN]]$AddChildNode(flat_tree[[i]])
  }
  return(assembled_tree)
}

to_data_frame <- function(xtree, active = NULL, values = c("nN", "obsN", "value", "error")) {
  # Returns data.frame, one row per node, 
  # with columns corresponding to requested
  # values. 
  # Arguments: 
    # xtree: full tree;
    # active: logical vector marking active nodes;
    # values: character vector specifying requested values;
  # Possible values are (defaults are marked with '*'):
    # "nN": node number, 
    # "parent_nN": parent node number,
    # "children_nN": children nodes,
    # "obsN": number of observation for a node, 
    # "value": predicted values for an observation that belongs to this node, 
    # "error": value of the error function for this node.
    # "split_conditions": split conditions to create this Node

  rv <- as.data.frame(lapply(values, function(v) sapply(xtree, function(x, value) { x[[v]] }, value = v)))
  colnames(rv) <- values
  rv
}

print_tree_old <- function(xtree, active = NULL){
  stopifnot(is.null(active) || length(active) == length(xtree))
  base::print(to_data_tree(xtree, active = active), "nN", "obsN", "value", "error")
}

plot_tree <- function(xtree, active = NULL){
  stopifnot(is.null(active) || length(active) == length(xtree))
  plot(to_data_tree(xtree, active = active))
}


tree_leaves <- function(xtree, active  = NULL, index = TRUE){
  # If isTRUE(index) returns index, i.e. nodes numbers, otherwise returns selector, i.e. logical vector of length(xtree)
  stopifnot(is.null(active) || length(active) == length(xtree))
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

tree_leaves_info <- function(xtree, active = NULL, as.dataframe = FALSE){
  leaves_index <- tree_leaves(xtree, active = active)
  leaves <- lapply(leaves_index, function(x) {
    xn <- xtree[[x]]
    xn[["defn"]] <- node_defn(xn, as_str = TRUE)
    xn <- xn[c("nN", "value", "error", "obsN", "defn")]
    xn
  })

  if(isTRUE(as.dataframe)){
    rv <- lapply(leaves, function(x) {
      # x <- within(x, split_conditions <- paste(split_conditions, collapse = " & "))
      x <- as.data.frame(x, stringsAsFactors = FALSE)
      row.names(x) <- NULL
      x
    })
    return(do.call(rbind, rv))
  }
  return(leaves)
}
    
predict_values <- function(xtree, X, active = NULL){
  leaves <- xtree[tree_leaves(xtree = xtree, active = active, index = TRUE)]
  leaves_values <- sapply(leaves, function(leaf) leaf$value)
  leaves_attrib <- lapply(leaves, function(leaf) eval(str2lang(paste(leaf[["split_conditions"]], collapse = " & ")), envir=X))
  leaves_attrib <- matrix(unlist(leaves_attrib), ncol = length(leaves), byrow = FALSE)
  Y <- leaves_attrib %*% leaves_values
  Y[, 1]
}


tree_info <- function(xtree, active = NULL) {

  # "tree_info" defines any nested subtree of the "xtree" 
  # and provides some "helper" information about this subtree
  # that allowes other functions to treet the subtree as "xtree"
  
  stopifnot(is.null(active) || length(active) == length(xtree))
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  active_nN <- seq_along(xtree)[active]
  leaves <- tree_leaves(xtree = xtree, active = active, index = TRUE)
  
  nd_sub_error <- cum_sub_error <- node_complexity <- rep(0, length(xtree))
  
  for (i in rev(active_nN)) {
    curr_node <- xtree[[i]]
    # nd_sub_error[i] <- with(curr_node, {if (i %in% leaves) 0 else error - sum(best_split$error)}) 
    # cum_sub_error[i] <- with(curr_node, {if (i %in% leaves) 0 else nd_sub_error[i] + (cum_sub_error[children_nN[1]] + cum_sub_error[children_nN[2]])})
    nd_sub_error[i] <- if (i %in% leaves) 0 else substitution_error(curr_node) 
    cum_sub_error[i] <- with(curr_node, {if (i %in% leaves) 0 else nd_sub_error[i] + (cum_sub_error[children_nN[1]] + cum_sub_error[children_nN[2]])})

    node_complexity[i] <- with(curr_node, {if (i %in% leaves) 1 else (node_complexity[children_nN[1]] + node_complexity[children_nN[2]])})
  }
  
  tmp <- rep(FALSE, length(xtree))
  tmp[leaves] <- TRUE
  return(list(active = active, node_substitution_error = nd_sub_error, node_complexity = node_complexity, cum_substitution_error = cum_sub_error, leaves = tmp))
}

path_list <- function(xtree, active = NULL){
  # Returns list of paths from all active leaves to the 
  # root node. Paths are expressed as vectors of node numbers
  # such that the first is always the root node, i.e. #1. Element
  # 1 in the list is a one-node path from and to the root node.
  # Element number 'n' will have a path from the node 'n' to the 
  # root node.
  # the last one is a leaf.
  # Paths are not padded, so vectors in the returned list have
  # different length.
  
  stopifnot(is.null(active) || length(active) == length(xtree))
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  # ti <- tree_info(xtree, active = active)
  
  path_L <- vector(mode = "list", length = sum(active))
  # max_path_length <- 1
  for (leaf in tree_leaves(xtree, active = active)) {
    path <- walk_up(xtree, leaf)
    path_L[[leaf]] <- rev(path)
    # if (length(path) > max_path_length) max_path_length <- length(path)
    path <- path[-1]
    for (i in seq_along(path)) {
      if (!is.null(path_L[[path[[i]]]])) {
        break
      } else {
        path_L[[path[[i]]]] <- rev(path[i:length(path)])
      }
    }
  }
  path_L
}

path_matrix <- function(xtree, active = NULL) {
  # Returns list from 'path_list' as a matrix
  # Order of rows is the same as in the xtree.
  # NB: this is not the same as ordering by 
  path_L <- path_list(xtree = xtree, active = active)
  max_path_length <- max(sapply(path_L, length))
  path_L <- lapply(path_L, function(x){ c(x, rep(NA, max_path_length - length(x))) })
  rv <- matrix(unlist(path_L), nrow = length(path_L), byrow = TRUE)
  rv
}


print_tree <- function(xtree, active = NULL){
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  tree_df <- to_data_frame(xtree)
  path_L <- path_list(xtree = xtree, active = active)
  max_path_length <- max(sapply(path_L, length))
  path_L <- lapply(path_L, function(x){ c(x, rep(NA, max_path_length - length(x))) })
  rv <- matrix(unlist(path_L), nrow = length(path_L), byrow = TRUE)
  in_order <- seq_len(length(xtree))
  out_order <- do.call(order, c(as.data.frame(rv), na.last = FALSE))
  rv <- rv[out_order, ]
  tree_df <- tree_df[out_order, ]
  active_out <- active[out_order]
  
  colwidth <- 3
  # stopifnot(colwidth %% 2 == 1)
  tab_length <- colwidth %/% 2 
  tab <- paste(rep(" ", tab_length), collapse = "")
  filler <- paste(rep(" ", colwidth), collapse = "")
  char_desc <- function(lev1, lev2 = NULL){
    if (is.null(lev2)) lev2 <- rep(NA, length(lev1))
    # browser()
    lev2_splits <- which(!is.na(lev2) & !duplicated(lev2))
    nsplits <- length(lev2_splits) / 2
    borders <- rep(filler, nrow(rv))
    for( i in seq_len(nsplits)) {
      s <- lev2_splits[2 * i - 1]
      e <- lev2_splits[2 * i]
      borders[s] <- paste(c(tab, "|", rep("-", tab_length)), collapse = "")
      borders[e] <- paste(c(tab, "*", rep("-", tab_length)), collapse = "")
      if ((e - s) > 1) borders[(s+1) : (e-1)] <- paste0(tab, "|", tab)
    }
    split_regions <- split(lev2_splits, factor(ceiling(seq_along(lev2_splits) / 2)))
    rv <- mapply(function(l, r, s){
      if (is.na(l)) return(tab)
      if (is.na(r)) return(with(xtree[[l]], split_conditions[length(split_conditions)]))
      return(s)
    }, lev1, lev2, borders)
    # browser()
    rv
  }
  level_L <- vector(mode = 'list', length = ncol(rv))
  for (i in seq_len(ncol(rv))) {
    lev1 <- rv[, i]
    if (i < ncol(rv)) lev2 <- rv[, i + 1] else lev2 <- NULL
    level_L[[i]] <- char_desc(lev1, lev2)
  }
  # browser()
  rv <- do.call(paste0, level_L)
  max_width <- max(sapply(rv, nchar))
  fmt <- paste0('%-', max_width, 's')
  rv <- sprintf(fmt, rv)
  rv <- cbind(Descr = rv, tree_df)
  return(rv[active_out, ])
}
