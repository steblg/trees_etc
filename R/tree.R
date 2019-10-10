# Inspired by:
# https://www.statworx.com/de/blog/coding-regression-trees-in-150-lines-of-code/

# Any split 's' devides split region into x <= s and x > s


rss_by_feature <- function(x, y, min_size = 2) {
  # browser()
  splits <- sort(unique(x))

  if(length(splits) == 1) return(NULL) # don't attempt to split if all x are the same
  rss_l <- lapply(splits, function(s){
    rss_f <- x <= s
    fn <- sum(rss_f) - 1
    nfn <- sum(!rss_f) - 1
    if(any(fn < min_size || nfn < min_size))
      rv <- c(NA, NA)
    else
      rv <- c(rss_l=var(y[rss_f]) * fn, rss_r=var(y[!rss_f]) * nfn)
    return(rv)
  })
  rss_v <- sapply(rss_l, sum)
  if(all(is.na(rss_v))) return(NULL) # not possible to split without breaking min_size requirement
  i <- which.min(rss_v)
  rv <- list(split = splits[i], rss = rss_l[[i]])
  return(rv)
}

gstree <- function(formula, input, model.control= list(minN = 20, minD = 5)) {
  # Do not split if number of observations is less or equal to minN
  # Do not split if drop in rss is less or equal to minD percent
  # First we implement minN only;

  input <- as.data.frame(input)

  formula_terms <- terms(formula, data = input)
  # stop if formula has no response
  stopifnot(attr(formula_terms, 'response') == 1)

  Y_label <- as.character(formula_terms)[2]
  Y <- model.response(model.frame(formula, input))

  # X_label <- attr(formula_terms, "terms.labels")
  X <- as.data.frame(model.matrix(formula, input)[, -1])
  X_label <- colnames(X)

  # # Status of any node can be either 'S' for 'Split', 'P' for 'Pass' and 'L' for a 'Leaf'
  # split_nodes <- list()

  # Tree to be grown
  rtree <- list()

  root_node <- list(
    nN = 1,                        # Node number
    status = 'S',
    parent = NULL,                 # Number of a parent Node
    split_conditions = "TRUE",       # List of split conditions to create this Node
    value = mean(Y),
    rss = var(Y) * (nrow(input) - 1),
    best_split = NULL,           # list with the following members:
    # split's feature and split's value, rss_l of left child node,
    # rss_r of the right child node
    obsN = nrow(input)           # Number of observations for this Node
  )

  rtree[[1]] <- root_node

  best_node_split <- function(node) {
    # Given the node number find the best possible split for this node
    # Returns a list with the following members:
    # "feature" as a number, "delta_rss" - change in rss due to the split,
    # "split" - value of the feature at which it is being split, "rss" - vector
    # with rss of the children nodes
    best_split <- node$best_split
    if(!is.null(best_split)) return(best_split)
    if(is.null(node$parent)){
      XX <- X
      YY <- Y
    } else {
      split_conditions_str <- paste(node$split_conditions, collapse = " & ")
      nodeF <- eval(str2lang(split_conditions_str), envir=X)  # not finished (*)
      # browser()
      XX <- X[nodeF, ]
      YY <- Y[nodeF]
    }
    # browser()
    # if( node$nN == 13) {
    #   possible_splits <- list()
    #   for(feature in seq_along(XX)){
    #     cat('Trying feature: ', X_label[feature], '\n')
    #     if(X_label[feature] == 'chas')
    #       debug(rss_by_feature)
    #     possible_splits[[feature]] <- rss_by_feature(XX[, feature], y = YY, min_size = model.control$minN)
    #   }
    # } else {
      possible_splits <- lapply(XX, rss_by_feature, y = YY, min_size = model.control$minN)
    # }
    tmp <- sapply(possible_splits, function(x) if(is.null(x)) NA else sum(x$rss))
    if(all(is.na(tmp))) return(NULL) # no allowed splits
    i <- which.min(tmp)
    return(c(list(feature = i), possible_splits[[i]]))
  }

  find_best_split <- function(nodes){
    split_effect <- sapply(nodes, function(i){
      return(rtree[[i]]$rss - sum(rtree[[i]]$best_split$rss))
    })
    return(nodes[which.max(split_effect)])
  }

  create_split <- function(i){
    node_to_split <- rtree[[i]]
    N <- length(rtree)

    children <- list()
    for(j in 1:2){
      child <- list(
        nN = N + j,
        parent = i
      )
      op <- c("<=", ">")
      # browser()
      split_condition_str <- with(node_to_split,
                                  c(split_conditions,
                                    paste(X_label[best_split$feature], op[j], best_split$split)
                                  )
      )
      child$split_conditions <- split_condition_str
      childF <- eval(str2lang(paste(split_condition_str, collapse = ' & ')), envir = X)
      child$obsN <- sum(childF)
      child$value <- mean(Y[childF])
      if(child$obsN <= model.control$minN) child$status <- 'L' else child$status <- 'S'

      # browser()
      child$rss <- with(node_to_split, best_split$rss[j])
      # browser()
      children[[j]] <- child
    }
    stopifnot((children[[1]]$obsN + children[[2]]$obsN) == node_to_split$obsN)
    return(children)
  }

  nodesN_to_split <- 1
  while(length(nodesN_to_split) > 0){

    cat("Nodes to split: ", paste(nodesN_to_split, collapse = ", "), "\n")
    # browser()
    split_info_not_available <- sapply(
      nodesN_to_split,
      function(i){is.null(rtree[[i]]$best_split)}
    )
    for(i in nodesN_to_split[split_info_not_available]){
      cat("Splitting node: ", i, "\n")
      # if(i == 17) browser()
      tmp <- best_node_split(rtree[[i]])
      if(is.null(tmp)) {
        rtree[[i]]$status <- 'L'
        nodesN_to_split <- seq_len(length(rtree))[sapply(rtree, function(x) x$status == 'S')]
      }
      rtree[[i]]$best_split <- tmp
    }

    if(length(nodesN_to_split) == 0 ) return(rtree)

    best_split <- find_best_split(nodesN_to_split)
    # if(is.null(best_split) || is.na(best_split)) browser()
    cat("Best split: ", paste(best_split, collapse = ", "), "\n")

    # update parent
    rtree[[best_split]]$status<- 'P'


    children <- create_split(best_split)
    rtree <- c(rtree, children)
    nodesN_to_split <- seq_len(length(rtree))[sapply(rtree, function(x) x$status == 'S')]

  }
  return(rtree)
}


# # Example:
# library(MASS)
# library(tree)
# data(Boston)
#
# set.seed(1)
# train <- sample(seq_len(nrow(Boston)), nrow(Boston)/2)
# foo <- gstree(medv ~ ., input = Boston[train, ])
