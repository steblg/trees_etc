
# Following recommendations from https://vctrs.r-lib.org/articles/s3-vector.html
# They are a bit excessive in our case, doing purely for learning purposes

# new_feature <- function(nm = character(), op = "<=") {
  # stopifnot(is.character(nm))
  # allowed_ops <- 
  # op <- match.arg(op, c("<=", "=="))
  # nop <- switch(

new_node <- function(n = 1L, parent = 0L, children = integer()){
  x <- list(
    nN = n,                         # Node number
    parent_nN = parent,             # Number of a parent Node
    children_nN = children,         # Children's numbers
    status = 'S'                    # Toggle, possible values 'S', for split-able, 'L' for a leaf node.
    # value = NULL,                   # Value assigned to the node (in regression trees, this is normally mean of node's observations)
    # error = NULL,                   # Value of the error function for this node
    # split_conditions = "TRUE",      # List of split conditions to create this Node
    # split_conditions_list = NULL,
    # best_split = NULL,              # list with the following members:
                                      # split's feature and split's value, rss_l of left child node,
                                      # error_right of the right child node
    # obsN = NULL                     # Number of observations for this Node
  )

  structure(
    x,
    class = 'te_node'
  )
}

validate_node <- function(x){
  with(x, {
    stopifnot(length(nN) > 0 && nN > 0)
    stopifnot(is.null(parent_nN) || parent_nN < nN)
    stopifnot(is.null(children_nN) || all(children_nN > nN))
  })
}

te_node <- function(n = integer(), parent = NULL, children = NULL){
  x <- new_node(n = n, parent = parent, children = children)
  validate_node(x)
  x
}

new_btree <- function(x = list()){ structure(x, class = "basic_tree") }

validate_tree <- function(x) { 
  stopifnot(all(sapply(x, class) == 'te_node')) 
  x
}

basic_tree <- function(x){
  return(validate_tree(new_btree(x)))
}

te_tree <- function(x = new_btree()){
  stopifnot(inherits(x, 'basic_tree'))
  ti <- tree_info(x)
  structure(
    ti,
    tree = x,
    class = 'te_tree'
  )
}
  
