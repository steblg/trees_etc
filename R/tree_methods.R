#To be reworked when using S3

node_defn <- function(node){
  slL <- list()
  for(e in node$split_conditions_list){
    slL[[e$feature]] <- c(slL[[e$feature]], list(e))
  }
  collapse_feature <- function(cL){ # elements of cL have the same 'feature'
      if(length(cL) == 1) return(cL)
      opF <- factor(sapply(cL, function(x) x$operation))
      L <- split(cL, opF)
      names(L) <- NULL
      collapse_F <- function(y) { # elements of y have same 'feature' and the same 'operation'
        if(y[[1]]$operation %in% c("==", '!=')) return(y)
        if(y[[1]]$operation == "<") return(list(feature = y[[1]]$feature, operation = y[[1]]$operation, split_value = min(sapply(y, function(x) x$split_value))))
        if(y[[1]]$operation == ">=") return(list(feature = y[[1]]$feature, operation = y[[1]]$operation, split_value = max(sapply(y, function(x) x$split_value))))
        stop(paste(paste("'", y[[1]]$operation, "'", sep =""), "is not supported"))
      }
      return(lapply(L, collapse_F))
  }
  return(lapply(slL, collapse_feature))
}

node_defn_str <- function(node){
  definition <- unlist(node_defn(node = node), recursive = FALSE, use.names = FALSE)
  return(paste(sapply(definition, function(x) with(x, paste("(", feature, operation, split_value, ")"))), collapse = " & "))
}
    
to_data_tree_node <- function(rtree_node, abbrev = TRUE) {
  library(data.tree)

  if(isTRUE(abbrev)){
    nm <- with(rtree_node, split_conditions[length(split_conditions)])
  } else {
    nm <- node_defn_str(rtree_node)
  }
  rv <- Node$new(nm)
  for (nm in names(rtree_node)) {
    rv[[nm]] <- rtree_node[[nm]]
  }
  rv
}

to_data_tree <- function(xtree, active = NULL){
  library(data.tree)
  
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  flat_tree <- lapply(xtree, to_data_tree_node)
  assembled_tree <- flat_tree[[1]]
  for(i in 2 : length(flat_tree)){
    if (active[i]) flat_tree[[flat_tree[[i]]$parent_nN]]$AddChildNode(flat_tree[[i]])
  }
  return(assembled_tree)
}

print_tree <- function(xtree){
  base::print(to_data_tree(xtree), "nN", "obsN", "value", "error")
}

plot_tree <- function(){
  plot(to_data_tree(xtree))
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

tree_leaves_info <- function(xtree, active = NULL, as.dataframe = FALSE){
  leaves_index <- tree_leaves(xtree, active = active)
  leaves <- lapply(leaves_index, function(x) {
    xn <- rtree_[[x]]
    xn[["defn"]] <- node_defn_str(xn)
    xn <- xn[c("nN", "split_conditions", "value", "error", "obsN", "defn")]
    xn
  })

  if(isTRUE(as.dataframe)){
    rv <- lapply(leaves_, function(x) {
      x <- within(x, split_conditions <- paste(split_conditions, collapse = " & "))
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

