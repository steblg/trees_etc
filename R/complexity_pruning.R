
pruning_sequence <- function(xtree, output = c("detailed", "summary"), pack = (output == "summary") ) {
  # Performs complexity pruning of the tree
  # Arguments:
  # xtree: basic tree to be pruned
  # For help: http://mlwiki.org/index.php/Cost-Complexity_Pruning
  
  info <- tree_info(xtree)
  info_list <- list(info)
  pruned_paths_list <- list(integer(0))
  alpha_list <- 0
  while ((an <- with(info, sum(active))) > 1) {
    # message("Number of active nodes: ", an)
    g_seq <- with(info, ifelse(active == TRUE, cum_substitution_error / (node_complexity - 1), NA))
    nnn <- order(g_seq, info$node_complexity, na.last = TRUE)[1]
    nn <- which.min(g_seq)
    stopifnot(nn == nnn)
    # message("Pruning: ", nn)
    # browser()
    prune_path <- c(tail(pruned_paths_list, 1)[[1]], nn)
    pruned_paths_list <- c(pruned_paths_list, list(prune_path))
    alpha <- g_seq[nn]
    info <- prune_node(xtree = xtree, node_n = nn, info = info)
    alpha_list <- c(alpha_list, alpha)
    info_list <- c(info_list, list(info))
  }
  output <- match.arg(output)
  if (output == "detailed") return(list(alpha = alpha_list, info = info_list))
  
  active_list <- lapply(info_list, function(x) x$active)
  if (pack == TRUE) {
    return(
      invisible(
        mapply(function(x, y, z) list(alpha = x, active = y, pruned = z), alpha_list, active_list, pruned_paths_list, SIMPLIFY = FALSE)
      )
    )
  } else {
    return(invisible(list(alpha = alpha_list, active = active_list, pruned = pruned_paths_list)))
  }
}

# cc_prune_tree <- function(xtree, alpha, prune_seq = NULL) {
  # # Returns a subtree, as a list of active nodes, 
  # # of the 'xtree' that corresponds to complexity multiplier 'alpha'
  # if (is.null(prune_seq)) prune_seq <- pruning_sequence(xtree = xtree, output = 'summary', pack = FALSE)
  # if (is.null(names(prune_seq))) {
    # alpha_list <- lapply(prune_seq, function(x) x$alpha)
    # active_list <- lapply(prune_seq, function(x) x$active)
    # prune_seq <- list(alpha = alpha_list, active = active_list)
  # }
  # alpha_list <- prune_seq$alpha
  # alpha <- sort(alpha)
  # ind <- findInterval(alpha, unlist(alpha_list), rightmost.closed = FALSE, all.inside = FALSE,left.open = FALSE)
  # return(prune_seq$active[ind])
# }

cc_prune_tree <- function(xtree, complexity_parameter) {
  # Returns a list of tree info, 
  # of the 'xtree' that corresponds to complexity multipliers 'complexity_parameter'
  prune_seq <- pruning_sequence(xtree = xtree, output = 'detailed')
  alpha_list <- prune_seq$alpha
  complexity_parameter <- sort(complexity_parameter)
  ind <- findInterval(complexity_parameter, unlist(alpha_list), rightmost.closed = FALSE, all.inside = FALSE,left.open = FALSE)
  return(prune_seq$info[ind])
}


# G. James, D. Witten, T. Hastie, R. Tibshirani, ISLR, 8th printing, 2017, page 309
# Algorithm 8.1 Building a regression tree

# 1. Use recursive binary splitting to grow a large tree on the training
# data, stopping only when each terminal node has fewer than some
# minimum number of observations.
# 2. Apply cost complexity pruning to the large tree in order to obtain a
# sequence of best subtrees, as a function of \alpha.
# 3. Use K-fold cross-validation to choose \αlpha. That is, divide the training
# observations into K folds. For each k = 1, . . . , K:
# (a) Repeat Steps 1 and 2 on all but the kth fold of the training data.
# (b) Evaluate the mean squared prediction error on the data in the
# left-out kth fold, as a function of \alpha.
# Average the results for each value of α, and pick α to minimize the
# average error.
# 4. Return the subtree from Step 2 that corresponds to the chosen value
# of \alpha.


# Or slightly more complete, and updated by me for clarity, explanation from SAS (https://documentation.sas.com/?docsetId=stathpug&docsetTarget=stathpug_hpsplit_details06.htm&docsetVersion=15.1&locale=en)
  # Finding the optimal subtree (or really complexity parameter) that does not overfit the training data is
  # performed by applying k-fold cross validation to cost-complexity pruning (Breiman et al. 1984; Zhang and
  # Singer 2010). The algorithm proceeds as follows after creating the sequence of subtrees and values by using
  # the entire set of training data :
 
  # Define cost complexity CC of a tree (T) as CC(T) = R(T) + \alpha * |T|, where R(T) represents T's error rate,
  # T| stands for a number of leaves on T while complexity parameter \alpha represents cost of each leaf. For 
  # classification trees misclassification rate is used as the error rate, R(T); for a continuous response
  # variable, the residual sum of squares (RSS), also called the sum of square errors (SSE), is used for the
  # error rate. Note that only the training data are used to evaluate cost complexity.  
  # Breiman et al. (1984) show that for each value of \alpha, there is a subtree of |T| that minimizes cost
  # complexity. When \alpha = 0, this is the full tree, T_0 . As \alpha increases, the corresponding subtree
  # becomes progressively smaller, and the subtrees are in fact nested. Then, at some value of \alpha, the root
  # node has the minimal cost complexity for any \alpha greater than or equal to that value. Because there are 
  # a finite number of possible subtrees, each subtree corresponds to an interval of values of ; that is,

  # [0, \alpha_1):          T_0
  # [\alpha_1, \alpha_2):   T_1
  
  # ....
  
  # [\alpha_n, Inf):        T_n

  # 1. Grow large tree on trainin data

  # 2. Apply cost complexity pruning to the large tree in order to obtain a sequence of best subtrees, 
  # as a function of \alpha. Define a sequence of \beta_i values as the geometric mean of the endpoints of the
  # [\alpha_i, \alpha_(i+1)) intervals to represent the intervals.

  # 3. Randomly divide the TRAINING observations into K approximately equal-sized parts, or folds.For each j-fold of the K folds, hold out the current fold for validation and use the remaining  k-1 folds for the training data in the following steps:

  # 3.1 Grow a large tree using k-1 folds data;

  # 3.2 Using the \beta_i values that are calculated in step 2, create a sequence of subtrees T_{j,i}, i.e. a subtree for each \beta_i, using cost-complexity pruning, but now using \beta_i as a fixed value \alpha for  and minimizing the cost-complexity of the tree, CC(T), to select a subtree T at each pruning step.

  # 3.3 For each \beta_i, set T (tree) to be the subtree T_{j,i} that has the minimum cost complexity from the sequence for the jth fold.

  # Calculate the error  for each  by using the current (jth) fold (the one omitted from the training).

  # Now the error rate can be averaged across folds, , and the  that has the smallest  is selected. The tree  from pruning the complete training data that corresponds to the selected  is the final selected subtree.

cvtune.new <- function(formula, input, k = 10, model.control = model_control()) { 
# remember that "input" here is a training data, not all available data!
# "cvtune" function is applied to the training data only!
# Returns

  input <- as.data.frame(input)

  formula_terms <- terms(formula, data = input)
  stopifnot(attr(formula_terms, 'response') == 1)

  Y <- model.response(model.frame(formula, input))

  X <- model.matrix(formula, input)
  ids <- attr(X, 'assign')
  ids <- which(ids > 0)
  X <- as.data.frame(X)[, ids, drop = FALSE]
  
  if (model.control$tree_type == 'regression') 
    test_err_func <- function(x, x_p){
      L <- sum(!is.na(x))
      return(sum((x - x_p)^2, na.rm = TRUE) / L)
      # var(x, na.rm=TRUE) * (L - 1)/L + (mean(x, na.rm = TRUE) - x_p)^2 
    }
  else
    test_err_func <- function(x, x_p) {
      L <- sum(!is.na(x))
      return(sum((x != x_p), na.rm = TRUE) / L)
    }
    
  
  # Following the explanation above:
  
  #1. 
  train_tree <- build_basic_tree(Y = Y, X = X, model.control = model.control)
  #2.
  train_prune_seq <- pruning_sequence(train_tree, output = "detailed")
  train_alpha <- train_prune_seq$alpha
  betas <- sqrt(train_alpha * c(train_alpha[-1], Inf))
  #3 
  fold_assignment <- sample(x = seq_len(k), size = nrow(input), replace = TRUE)
  process_fold <- function(test_fold, Y, X, model.control){
    train_fold <- !test_fold
    #3.1
    fold_tree <- build_basic_tree(Y = Y[train_fold], X = X[train_fold, , drop = FALSE], model.control = model.control)
    #3.2
    fold_prune_seq <- cc_prune_tree(xtree = fold_tree, complexity_parameter = betas)
    #3.3
    test_err_list <- sapply(fold_prune_seq, function(prune_info, tree, x_test, y_test) {
        y_pred <- predict_values(xtree = tree, X = x_test, active = prune_info$active)
        test_err <- test_err_func(y_test, y_pred)
        return(test_err)
      }, tree = fold_tree, x_test = X[test_fold, , drop = FALSE], y_test = Y[test_fold]
    )
    return(test_err_list)
  }

  alpha_ranges_num <- length(train_alpha)
  betas_num <- length(betas)
  stopifnot(alpha_ranges_num == betas_num)
  fold_err_vec <- vector(mode = "numeric", length = k * alpha_ranges_num)
  for (i in 1:k) {
    test_fold <- fold_assignment == i
    fold_err <- process_fold(test_fold = test_fold, Y = Y, X = X, model.control = model.control)
    fs <- (i - 1) * betas_num + 1
    fe <- i * betas_num
    fold_err_vec[fs : fe] <- fold_err
  }
  # calculate average error per train_alpha range
  betas_ind <- rep(seq_len(betas_num), k)
  # browser()
  avg_error <- aggregate(fold_err_vec, by = list(betas = betas_ind), FUN = function(...) mean(..., na.rm=TRUE))[, 2]
  return(list(tree = train_tree, prune_seq = train_prune_seq, betas = betas, avg_error = avg_error))
}


cvtune.old <- function(formula, input, k = 10, model.control = model_control()) { 
# remember that "input" here is a training data, not all available data!
# "cvtune" function is applied to the training data only!
# Returns

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

cvtune <- cvtune.new
  
