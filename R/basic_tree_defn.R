
#var_error <- function(x) var(x)

getmode <- function(x) {
  uv <- unique(x)
  return(uv[which.max(table(uv))])
}

gini_impurity <- function(x) {
  probs <- table(x) / length(x)
  return(1 - sum(probs^2))
}

entropy <- function(x) {
  probs <- table(x) / length(x)
  return(sum(-probs * log2(probs)))
}

misclass <- function(x, x_p) {
  L <- sum(!is.na(x))
  stopifnot(L > 0)
  if (missing(x_p)) {
    probs <- table(x) / L
    return(1 - max(probs))
  } else {
    return(sum(x != x_p, na.rm = TRUE) / L)
  }
}

feature_split_num <- function(x){
  # This is not the most efficient, or best, way of defining split points.
  # ToDo: come up with a better one

  breakpoints <- sort(unique(x), na.last = NA)
  breakpoints <- (breakpoints[-length(breakpoints)] + breakpoints[-1]) / 2 # select midpoints between observations
  op <- '<'
  return(list(s = breakpoints, op=op))
}

feature_split_cat <- function(x){
  xl <- levels(droplevels(x))
  unique_sets <- sample_space(xl, complimentary = FALSE)
  op <- '%in%'
  return(list(s=unique_sets, op=op))
}

feature_split <- function(x){
  # Create list of possible splits
  return(if (is.factor(x)) feature_split_cat(x) else feature_split_num(x))
}

split_predictor_region <- function(x, y, err_FUN, min_size) {
  # For every possible split along predictor 'x'
  # calculate "errors" of the resulting nodes;

  splits <- feature_split(x)
  op <- splits$op
  if(length(splits) == 1) return(NULL) # don't attempt to split if all x are the same
  N <- length(y)
  split_list <- lapply(splits$s, function(s){
    split_filter <- do.call(op, list(x, s))
    obs_N <- table(split_filter)
    if(any(obs_N[2] < min_size || obs_N[1] < min_size))
      rv <- list(err = c(NA, NA), probs = rev(obs_N) / N)
    else
      rv <- list(err = c(err_FUN(y[split_filter]), err_FUN(y[!split_filter])), probs = rev(obs_N) / N)
    # browser()
    return(rv)
  })
  # browser()
  error_value <- sapply(split_list, function(X)sum(X$err * X$probs))
  if(all(is.na(error_value))) return(NULL) # not possible to split without breaking min_size requirement
  i <- which.min(error_value)
  rv <- list(split = splits$s[i], op = op, error = split_list[[i]])
  return(rv)
}

bestNodeSplit <- function(node, X, Y, model.control) {
  # Given the node find the best possible split for this node
  # 'X' is model matrix
  # 'Y' is model response
  best_split <- node$best_split
  if(!is.null(best_split)) return(best_split)
  X_label <- colnames(X)
  if(is.null(node$parent_nN)){
    XX <- X
    YY <- Y
  } else {
    split_conditions_str <- paste(node$split_conditions, collapse = " & ")
    nodeF <- eval(str2lang(split_conditions_str), envir=X)  # not finished (*)
    XX <- X[nodeF, ]
    YY <- Y[nodeF]
  }
  possible_splits <- lapply(XX, split_predictor_region, y = YY, err_FUN = model.control$error_func, min_size = model.control$min_node_size)
  # browser()
  tmp <- sapply(possible_splits, function(x) if(is.null(x)) NA else sum(x$error$err * x$error$probs))
  if(all(is.na(tmp))) return(NULL) # no allowed splits
  i <- which.min(tmp)
  return(c(list(feature = X_label[i]), possible_splits[[i]]))
}

substitution_error <- function(node){
  bse <- node$best_split$error
  return(node$obsN * (node$error - sum(bse$err * bse$probs)))
}

findBestNodeToSplit <- function(nodes, xtree){
  err_delta <- sapply(nodes, function(i){
    # # return(xtree[[i]]$error - sum(xtree[[i]]$best_split$error))
    # return(with(xtree[[i]], obsN * (error - sum(with(best_split$error, err * probs)))))
    return(substitution_error(xtree[[i]]))
  })
  return(nodes[which.max(err_delta)])
}

createSplit <- function(node_to_split, tree_length, X, Y, model.control){
  # 'X' is model matrix
  # 'Y' is model response

  # node_to_split <- xtree[[i]]
  N <- tree_length

  children <- list()
  bs <- node_to_split$best_split
  bsop <- bs$op
  if(bsop == "<") {
    nop <- ">="
  } else if(bsop == "==") {
    nop <- "!="
  } else {
    stop(paste("Operation", paste("'", bsop, "'", sep = ""), "not supported"))
  }
  op <- c(bsop, nop)
  for(j in 1:2){
    child <- list(
      nN = N + j,
      parent_nN = node_to_split$nN
    )
    add_condition <- list(feature = bs$feature, operation = op[j], split_value = bs$split)
    split_condition_str <- c(node_to_split$split_conditions, do.call(paste, (add_condition)))
    childF <- eval(str2lang(paste(split_condition_str, collapse = ' & ')), envir = X)

    child <- within(child, {
      split_conditions      <- split_condition_str
      split_conditions_list <- c(node_to_split$split_conditions_list, list(add_condition))
      obsN                  <- sum(childF)
      error                 <- node_to_split$best_split$error$err[j]
      value                 <- model.control$value_func(Y[childF])
      status                <- 'S'
    })
    child$best_split <- tmp <- bestNodeSplit(child, X = X, Y = Y, model.control = model.control)
    if (is.null(tmp)) child <- markNodeLeaf(child)
    children[[j]] <- child
  }
  stopifnot((children[[1]]$obsN + children[[2]]$obsN) == node_to_split$obsN)
  return(children)
}

markNodeLeaf <- function(node){
  node$status <- 'L'
  node['best_split'] <- list(NULL)
  node['children_nN'] <- list(NULL)
  return(node)
}


build_basic_tree <- function(Y, X, model.control = model_control()) {
  # 'X' is model matrix
  # 'Y' is model response
  
  errorFUN <- model.control$error_func
  min_size <- model.control$min_node_size
  X_label <- colnames(X)

  # Tree to be grown
  rtree_ <- list()

  root_node <- list(
    nN = 1,                        # Node number
    status = 'S',
    parent_nN = NULL,                 # Number of a parent Node
    split_conditions = "TRUE",
    split_conditions_list = NULL,     # List of split conditions to create this Node
    value = mean(Y),
    error = errorFUN(Y),
    best_split = NULL,           # list with the following members:
    # split's feature and split's value, rss_l of left child node,
    # error_right of the right child node
    obsN = nrow(X)           # Number of observations for this Node
  )

  tmp <- bestNodeSplit(root_node, X = X, Y = Y, model.control = model.control)
  if(is.null(tmp)) { # no allowed splits
    stop("Tree can't be build under current conditions, try changing model.control argument")
  }
  root_node['best_split'] <- list(tmp)
  rtree_[[1]] <- root_node

  # error value of the entire_tree
  errorVal <- with(root_node, error * obsN)

  nodesN_to_split <- 1
  while(length(nodesN_to_split) > 0){

    # cat("Nodes that can be split: ", paste(nodesN_to_split, collapse = ", "), "\n")

    best_split <- findBestNodeToSplit(nodesN_to_split, xtree = rtree_)
    # cat("Best split: ", paste(best_split, collapse = ", "), "\n")
    # update parent
    nodeToSplit <- rtree_[[best_split]]
    # deltaError <- with(nodeToSplit, error - sum(best_split$error))
    deltaError <- substitution_error(nodeToSplit)
    if(100 * deltaError / errorVal[length(errorVal)] < model.control$min_sub_err){
      for(i in nodesN_to_split) rtree_[[i]] <- markNodeLeaf(rtree_[[i]])
      break
    } else {
      kids <- createSplit(node_to_split = nodeToSplit, tree_length = length(rtree_), X = X, Y = Y, model.control)
      rtree_[[best_split]] <- within(nodeToSplit, {status <- 'P'; children_nN <- sapply(kids, function(x) x$nN)})
      rtree_ <- c(rtree_, kids)
      errorVal <- errorVal - nodeToSplit$error * nodeToSplit$obsN + sum(sapply(kids, function(x) x$error * x$obsN))
    }
    nodesN_to_split <- seq_len(length(rtree_))[sapply(rtree_, function(x) x$status == 'S')]
  }
  return(rtree_)
}
