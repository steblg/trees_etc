rss <- function(x) var(x) * (length(x) - 1)
# rss <- function(x) var(x) * length(x)

basic_tree <- function(formula, input, model.control= model_control()) {
    input <- as.data.frame(input)

    formula_terms <- terms(formula, data = input)
    # stop if formula has no response
    stopifnot(attr(formula_terms, 'response') == 1)

#    Y_label <- as.character(formula_terms)[2]
    Y <- model.response(model.frame(formula, input))

    # X_label <- attr(formula_terms, "terms.labels")
    X <- model.matrix(formula, input)
    ids <- attr(X, 'assign')
    ids <- which(ids > 0)
    X <- as.data.frame(X)[, ids, drop = FALSE]
    return(build_basic_tree(Y = Y, X = X, model.control = model.control))
}

model_control <- function(minS = 20, minD = 5, error = c("deviance", "gini")) {
  # Do not split if number of observations is less or equal to minS
  # Do not split if drop in rss is less or equal to minD percent
  error <- match.arg(error)
  errorFun <- switch(error, deviance = rss, gini = function(...) stop("Not implemented"))
  return(list(minS = minS, minD = minD, errorFun = errorFun))
}

build_basic_tree <- function(Y, X, model.control = model_control()) {

  if( is.factor(Y) ){
    stopifnot( ! is.factor(Y) ) # Classification trees aren't implemented yet.
  } else {
    minFUN <- model.control[['errorFun']]
  }
  min_size <- model.control$minS
  X_label <- colnames(X)

  split_along_predictor <- function(x, y) {
    if(is.factor(x)){
      # stopifnot(require("arrangements"))  # Combinatorics support package
      splits <- levels(x) # this must be changed, too simplistic and not complete
      op <- '%in%'
    } else {
      splits <- sort(unique(x))
      splits <- (splits[-length(splits)] + splits[-1]) / 2 # select midpoints between observations
      op <- '<'
    }
    if(length(splits) == 1) return(NULL) # don't attempt to split if all x are the same
    split_list <- lapply(splits, function(s){
      split_filter <- do.call(op, list(x, s))
      filter_in <- sum(split_filter) - 1
      filter_out <- sum(!split_filter) - 1
      if(any(filter_in < min_size || filter_out < min_size))
        rv <- c(NA, NA)
      else
        rv <- c(minFUN(y[split_filter]), minFUN(y[!split_filter]))
      return(rv)
    })
    error_value <- sapply(split_list, sum)
    if(all(is.na(error_value))) return(NULL) # not possible to split without breaking min_size requirement
    i <- which.min(error_value)
    rv <- list(split = splits[i], op = op, error = split_list[[i]])
    return(rv)
  }

  bestNodeSplit <- function(node) {
    # Given the node find the best possible split for this node
    # Returns a list with the following members:
    # "feature" as a number, "delta_rss" - change in rss due to the split,
    # "split" - value of the feature at which it is being split, "rss" - vector
    # with rss of the children nodes
    best_split <- node$best_split
    if(!is.null(best_split)) return(best_split)
    if(is.null(node$parent_nN)){
      XX <- X
      YY <- Y
    } else {
      split_conditions_str <- paste(node$split_conditions, collapse = " & ")
      nodeF <- eval(str2lang(split_conditions_str), envir=X)  # not finished (*)
      XX <- X[nodeF, ]
      YY <- Y[nodeF]
    }
    possible_splits <- lapply(XX, split_along_predictor, y = YY)
    tmp <- sapply(possible_splits, function(x) if(is.null(x)) NA else sum(x$error))
    if(all(is.na(tmp))) return(NULL) # no allowed splits
    i <- which.min(tmp)
    return(c(list(feature = i), possible_splits[[i]]))
  }

  findBestNodeToSplit <- function(nodes){
    split_effect <- sapply(nodes, function(i){
      return(rtree_[[i]]$error - sum(rtree_[[i]]$best_split$error))
    })
    return(nodes[which.max(split_effect)])
  }

  markNodeLeaf <- function(node){
    node$status <- 'L'
    node['best_split'] <- list(NULL)
    node['children_nN'] <- list(NULL)
    return(node)
  }

  createSplit <- function(i){
    node_to_split <- rtree_[[i]]
    N <- length(rtree_)

    children <- list()
    op <- with(node_to_split, best_split$op)
    if(op == "<") {
      nop <- ">="
    } else if(op == "==") {
      nop <- "!="
    } else {
      stop(paste("Operation", paste("'", op, "'", sep = ""), "not supported"))
    }
    op <- c(op, nop)
    for(j in 1:2){
      child <- list(
        nN = N + j,
        parent_nN = i
      )
      add_condition <- with(node_to_split, list(feature=X_label[best_split$feature], operation = op[j], split_value=best_split$split))
      split_condition_str <- c(node_to_split$split_conditions, with(add_condition, paste(feature, operation, split_value)))
      childF <- eval(str2lang(paste(split_condition_str, collapse = ' & ')), envir = X)

      child <- within(child, {
        split_conditions      <- split_condition_str
        split_conditions_list <- c(node_to_split$split_conditions_list, list(add_condition))
        obsN                  <- sum(childF)
        error                 <- node_to_split$best_split$error[j]
        value                 <- mean(Y[childF])
        status                <- 'S'
      })
      child$best_split <- tmp <- bestNodeSplit(child)
      if(is.null(tmp)) child <- markNodeLeaf(child)
      children[[j]] <- child
    }
    stopifnot((children[[1]]$obsN + children[[2]]$obsN) == node_to_split$obsN)
    return(children)
  }

  # Tree to be grown
  rtree_ <- list()
  # error value of the model
  errorVal <- minFUN(Y)



  root_node <- list(
    nN = 1,                        # Node number
    status = 'S',
    parent_nN = NULL,                 # Number of a parent Node
    split_conditions = "TRUE",       
    split_conditions_list = NULL,     # List of split conditions to create this Node
    value = mean(Y),
    error = minFUN(Y),
    best_split = NULL,           # list with the following members:
    # split's feature and split's value, rss_l of left child node,
    # error_right of the right child node
    obsN = nrow(X)           # Number of observations for this Node
  )

  tmp <- bestNodeSplit(root_node)
  if(is.null(tmp)) { # no allowed splits
    stop("Tree can't be build under current conditions, try changing model.control argument")
  }
  root_node['best_split'] <- list(tmp)
  rtree_[[1]] <- root_node

  nodesN_to_split <- 1
  while(length(nodesN_to_split) > 0){

    # cat("Nodes that can be split: ", paste(nodesN_to_split, collapse = ", "), "\n")

    best_split <- findBestNodeToSplit(nodesN_to_split)
    # cat("Best split: ", paste(best_split, collapse = ", "), "\n")
    # update parent
    nodeToSplit <- rtree_[[best_split]]
    deltaError <- with(nodeToSplit, error - sum(best_split$error))
    if(100 * deltaError / errorVal[length(errorVal)] < model.control$minD){
      for(i in nodesN_to_split) rtree_[[i]] <- markNodeLeaf(rtree_[[i]])
      break
    } else {
      kids <- createSplit(best_split)
      rtree_[[best_split]] <- within(nodeToSplit, {status <- 'P'; children_nN <- sapply(kids, function(x) x$nN)})
      rtree_ <- c(rtree_, kids)
      errorVal <- c(errorVal, sum(sapply(rtree_, function(x){
        if(x$status == 'P') return(0) else return(x$error)
      })))
    }
    nodesN_to_split <- seq_len(length(rtree_))[sapply(rtree_, function(x) x$status == 'S')]
  }

  # rv <- list(
    # tree = rtree_,
    # errVal = errorVal
  # )
  return(rtree_)
}
