#To be reworked when using S3

nodeDefn <- function(node){
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

nodeDefnStr <- function(node){
  definition <- unlist(nodeDefn(node = node), recursive = FALSE, use.names = FALSE)
  return(paste(sapply(definition, function(x) with(x, paste("(", feature, operation, split_value, ")"))), collapse = " & "))
}
    
toDataTreeNode <- function(rtree_node, abbrev = TRUE) {
  library(data.tree)

  if(isTRUE(abbrev)){
    nm <- with(rtree_node, split_conditions[length(split_conditions)])
  } else {
    nm <- nodeDefnStr(rtree_node)
  }
  rv <- Node$new(nm)
  for (nm in names(rtree_node)) {
    rv[[nm]] <- rtree_node[[nm]]
  }
  rv
}

toDataTree <- function(xtree, active = NULL){
  library(data.tree)
  
  if (is.null(active)) active <- rep(TRUE, length(xtree))
  flat_tree <- lapply(xtree, toDataTreeNode)
  assembled_tree <- flat_tree[[1]]
  for(i in 2 : length(flat_tree)){
    if (active[i]) flat_tree[[flat_tree[[i]]$parent_nN]]$AddChildNode(flat_tree[[i]])
  }
  return(assembled_tree)
}

printTree <- function(xtree){
  base::print(toDataTree(xtree), "nN", "obsN", "value", "error")
}

plotTree <- function(){
  plot(toDataTree(xtree))
}


treeLeaves <- function(xtree, active  = NULL, index = TRUE){
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

printLeaves <- function(xtree, active = NULL, as.dataframe = FALSE){
  leaves_index <- treeLeaves(xtree, active = active)
  leaves <- lapply(leaves_index, function(x) {
    xn <- rtree_[[x]]
    xn[["defn"]] <- nodeDefnStr(xn)
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
    
predictValues <- function(xtree, X, active = NULL){
  leaves <- xtree[treeLeaves(xtree = xtree, active = active, index = TRUE)]
  leaves_values <- sapply(leaves, function(leaf) leaf$value)
  leaves_attrib <- lapply(leaves, function(leaf) eval(str2lang(leaf[["split_conditions"]]), envir=X))
  leaves_attrib <- matrix(unlist(leaves_attrib), ncol = length(leaves), byrow = FALSE)
  Y <- leaves_attrib %*% leaves_values
  Y
}
