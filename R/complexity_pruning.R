
#basic_tree <- function(formula, input, model.control= list(minS = 20, minD = 5, error = c("deviance", "gini"))) 

pruning_sequence <- function(xtree) {
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
    message("Number of active nodes: ", an)
    g_seq <- with(info, ifelse(active == TRUE, substitution_error / (node_complexity - 1), NA))
    nn <- which.min(g_seq)
    message("Pruning: ", nn)
    alpha <- g_seq[nn]
    info <- prune_node(xtree = xtree, node_n = nn, info = info)
    alpha_list <- c(alpha_list, alpha)
    info_list <- c(info_list, list(info))
  }
  alpha_list <- c(sqrt(alpha_list[-length(alpha_list)] * alpha_list[-1]), alpha_list[length(alpha_list)])
  complexity_penalty_list <- alpha_list * sapply(info_list, function(x) x$complexity[1])
  active_list <- lapply(info_list, function(x) x$active)
  return(
    invisible(
      mapply(function(x, y, z) list(alpha = x, cp = y, active = z), alpha_list, complexity_penalty_list, active_list, SIMPLIFY = FALSE)
    )
  )
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
# remember that "input" here is a training data, not all available !
# "cvtune" function is applied to the training data only!

  input <- as.data.frame(input)

  formula_terms <- terms(formula, data = input)
  stopifnot(attr(formula_terms, 'response') == 1)

  Y <- model.response(model.frame(formula, input))

  X <- model.matrix(formula, input)
  ids <- attr(X, 'assign')
  X <- as.data.frame(X[, ids, drop = FALSE])

  fold_assignment <- sample(x = 1:k, size = nrow(input), replace = TRUE)
  folds <- split(seq_len(nrow(input), f = fold_assignment))
  L <- lapply(folds, function(test_fold, Y, X, model.control){
      train_fold <- seq_len(nrow(X))[-test_fold]
      fold_tree <- build_basic_tree(Y = Y[train_fold], X = X[traind_fold, , drop = FALSE], model.control = model.control)
      fold_prune_seq <- pruning_sequence(fold_tree)
      test_err_list <- lapply(fold_prune_seq, function(prune_info, tree, x_test, y_test) {
          y_pred <- predict_values(xtree = tree, X = x_test, active = prune_info$active)
          test_err <- (y_pred - y_test)^2 / length(y_test)
        }, tree = fold_tree, x_test = X[test_fold, , drop = FALSE], y_test = Y[test_fold]
      )
    },
    Y = Y, X = X, model.control = model.control
  )