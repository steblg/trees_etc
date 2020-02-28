
model_control <- function(tree_type = NULL, min_node_size = NULL, min_sub_err = NULL) {
  # tree_type can be either NULL, "regression" or "classification"
  # Do not split if number of observations is less or equal to min_node_size
  # Do not split if substitution error is less or equal to min_sub_err percent
  if (is.null(tree_type)) tree_type <- "regression"
  stopifnot(tree_type %in% c("regression", "classification"))
  default_controls <- list(
    tree_type = if (!is.null(tree_type)) tree_type else "regression",
    min_node_size = if (!is.null(min_node_size)) min_node_size else 10,
    min_sub_err = if (!is.null(min_sub_err)) min_sub_err else 5
  )
  if (tree_type == 'regression')
    func_list <- list(
      value_func = function(x) mean(x, na.rm = TRUE),
      error_func = function(x) {
        L <- sum(!is.na(x))
        return((1 - 1/L) * var(x, na.rm = TRUE))
      },
      prune_error_func = function(x) var(x, na.rm = TRUE)
    )
  else
    func_list <- list(
      value_func = getmode,
      error_func = gini_impurity,
      prune_error_func = misclass
    )
    
  controls <- c(default_controls, func_list)
  return(controls)
}


basic_tree <- function(formula, input, model.control= model_control()) {
    input <- as.data.frame(input)

    formula_terms <- terms(formula, data = input)
    # stop if formula has no response
    stopifnot(attr(formula_terms, 'response') == 1)

#    Y_label <- as.character(formula_terms)[2]
    Y <- model.response(model.frame(formula, input))
    
    if (any(is.na(Y))) stop("Response variable must not have NA values")
    
    if (is.factor(Y)) 
      stopifnot(model.control$tree_type == 'classification')
    else
      stopifnot(model.control$tree_type == 'regression')

    # X_label <- attr(formula_terms, "terms.labels")
    X <- model.matrix(formula, input)
    ids <- attr(X, 'assign')
    ids <- which(ids > 0)
    X <- as.data.frame(X)[, ids, drop = FALSE]
    return(build_basic_tree(Y = Y, X = X, model.control = model.control))
}
