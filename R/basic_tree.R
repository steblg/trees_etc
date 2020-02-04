
model_control <- function(...) {
  # Do not split if number of observations is less or equal to min_node_size
  # Do not split if substitution error is less or equal to min_sub_err percent
  default_controls <- list(
    min_node_size = 10,
    min_sub_err = 5,
    error_func = 'variance'
  )
  requested_controls <- list(...)
  controls <- default_controls
  if (length(requested_controls) > 0) {
    for(nm in names(requested_controls)) 
      if (nm %in% names(controls)) 
        controls[[nm]] <- requested_controls[[nm]]
  }
  
  controls <- within( controls, errorFun <- switch(error_func, variance = var_error, gini = gini_impurity, entropy = entropy))
  return(controls)
}

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
