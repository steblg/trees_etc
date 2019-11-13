
#basic_tree <- function(formula, input, model.control= list(minS = 20, minD = 5, error = c("deviance", "gini"))) 

pruning_sequence <- function(xtree, pack = TRUE) {
  # Performs complexity pruning of the tree
  # Arguments:
  # xtree: basic tree to be pruned
  # Returns list of elements, one per prune step:
  # alpha: regularization parameter (complexity multiplier);
  # cp: complexity_penalty, i.e. alpha * number of leaves in a pruned tree
  # active: active nodes in the pruned tree
  
  info <- tree_info(xtree)
  info_list <- list(info)
  alpha_list <- 0
  while ((an <- sum(info$active)) > 1) {
    # message("Number of active nodes: ", an)
    g_seq <- with(info, ifelse(active == TRUE, substitution_error / (node_complexity - 1), NA))
    nnn <- order(g_seq, info$node_complexity, na.last = TRUE)[1]
    nn <- which.min(g_seq)
    stopifnot(nn == nnn)
    # message("Pruning: ", nn)
    alpha <- g_seq[nn]
    info <- prune_node(xtree = xtree, node_n = nn, info = info)
    alpha_list <- c(alpha_list, alpha)
    info_list <- c(info_list, list(info))
  }
#  alpha_list <- c(sqrt(alpha_list[-length(alpha_list)] * alpha_list[-1]), alpha_list[length(alpha_list)])

#  complexity_penalty_list <- alpha_list * sapply(info_list, function(x) x$complexity[1])
  active_list <- lapply(info_list, function(x) x$active)
  if (pack == TRUE) {
    return(
      invisible(
#        mapply(function(x, y, z) list(alpha = x, cp = y, active = z), alpha_list, complexity_penalty_list, active_list, SIMPLIFY = FALSE)
        mapply(function(x, y) list(alpha = x, active = y), alpha_list, active_list, SIMPLIFY = FALSE)
      )
    )
  } else {
    # return(invisible(list(alpha = alpha_list, cp = complexity_penalty_list, active = active_list)))
    return(invisible(list(alpha = alpha_list, active = active_list)))
  }
}

cc_prune_tree <- function(xtree, alpha, prune_seq = NULL) {
  if (is.null(prune_seq)) prune_seq <- pruning_sequence(xtree = xtree, pack = FALSE)
  if (is.null(names(prune_seq))) {
    alpha_list <- lapply(prune_seq, function(x) x$alpha)
    active_list <- lapply(prune_seq, function(x) x$active)
    prune_seq <- list(alpha = alpha_list, active =active_list)
  }
  alpha_list <- prune_seq$alpha
  alpha <- sort(alpha)
  ind <- findInterval(alpha, unlist(alpha_list), rightmost.closed = FALSE, all.inside = FALSE,left.open = FALSE)
  return(prune_seq$active[ind])
}

# G. James, D. Witten, T. Hastie, R. Tibshirani, ISLR, 8th printing, 2017, page 309
# Algorithm 8.1 Building a regression tree

# 1. Use recursive binary splitting to grow a large tree on the training
# data, stopping only when each terminal node has fewer than some
# minimum number of observations.
# 2. Apply cost complexity pruning to the large tree in order to obtain a
# sequence of best subtrees, as a function of α.
# 3. Use K-fold cross-validation to choose α. That is, divide the training
# observations into K folds. For each k = 1, . . . , K:
# (a) Repeat Steps 1 and 2 on all but the kth fold of the training data.
# (b) Evaluate the mean squared prediction error on the data in the
# left-out kth fold, as a function of α.
# Average the results for each value of α, and pick α to minimize the
# average error.
# 4. Return the subtree from Step 2 that corresponds to the chosen value
# of α.

cvtune <- function(formula, input, k = 10, model.control) { 
# remember that "input" here is a training data, not all available data!
# "cvtune" function is applied to the training data only!

  input <- as.data.frame(input)

  formula_terms <- terms(formula, data = input)
  stopifnot(attr(formula_terms, 'response') == 1)

  Y <- model.response(model.frame(formula, input))

  X <- model.matrix(formula, input)
  ids <- attr(X, 'assign')
  ids <- which(ids > 0)
  X <- as.data.frame(X)[, ids, drop = FALSE]
  
  # Following the algorithm above:
  
  #1. 
  train_tree <- build_basic_tree(Y = Y, X = X, model.control = model.control)
  #2.
  train_prune_seq <- pruning_sequence(train_tree, pack = FALSE)
  train_alpha <- train_prune_seq$alpha
  train_alpha_vals <- c(sqrt(train_alpha[-length(train_alpha)] * train_alpha[-1]), train_alpha[length(train_alpha)])
  #3
  fold_assignment <- sample(x = 1:k, size = nrow(input), replace = TRUE)
  # folds <- split(seq_len(nrow(input), f = fold_assignment))
  # For each training fold calc pruning sequence and a series of corresponding prediction errors
  #L <- lapply(folds, function(test_fold, Y, X, model.control){
  process_fold <- function(test_fold, Y, X, model.control){
    train_fold <- !test_fold
    fold_tree <- build_basic_tree(Y = Y[train_fold], X = X[train_fold, , drop = FALSE], model.control = model.control)
    fold_prune_seq <- pruning_sequence(fold_tree, pack = TRUE)

    test_err_list <- lapply(fold_prune_seq, function(prune_info, tree, x_test, y_test) {
        y_pred <- predict_values(xtree = tree, X = x_test, active = prune_info$active)
        test_err <- sum((y_pred - y_test)^2) / (length(y_test) - 1)
        return(list(alpha = prune_info$alpha, err = test_err))
      }, tree = fold_tree, x_test = X[test_fold, , drop = FALSE], y_test = Y[test_fold]
    )
    return(test_err_list)
  }
  alpha_ranges_num <- length(train_alpha)
  fold_err_vec <- vector(mode = "numeric", length = k * alpha_ranges_num)
  for (i in 1:k) {
    test_fold <- fold_assignment == i
    fold_info <- process_fold(test_fold = test_fold, Y = Y, X = X, model.control = model.control)
    # find how alphas of this fold pruning sequence fit into ranges of train_alpha
    alpha_assignment <- findInterval(
      sapply(fold_info, function(x) x$alpha), 
      train_alpha, 
      rightmost.closed = FALSE, all.inside = FALSE,left.open = FALSE
    )
    # calculate this fold's average prediction error for each train_alpha range
    fold_err <- rep(NA_real_, length(train_alpha))
    agg_errs <- aggregate(sapply(fold_info, function(x) x$err), by = list(alpha_assignment), FUN = mean)
    fold_err[agg_errs[, 1]] <- agg_errs[, 2]
    fs <- (i - 1) * alpha_ranges_num + 1
    fe <- i * alpha_ranges_num
    fold_err_vec[fs : fe] <- fold_err
  }
  # calculate average error per train_alpha range
  alpha_range_ind <- rep(1:alpha_ranges_num, k)
  # browser()
  avg_error <- aggregate(fold_err_vec, by = list(alpha_range = alpha_range_ind), FUN = function(...) mean(..., na.rm=TRUE))[, 2]
  # browser()
  return(list(tree = train_tree, prune_seq = train_prune_seq, alpha = train_alpha, avg_error = avg_error))
}

  
