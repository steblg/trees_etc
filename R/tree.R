

# ToDo:
#   1. Handling of NA's. For now one has to deal with them (na.omit) prior to using the method;
#   2. Feature importance
#   3. Complexity pruning
#   4. Cross-validation


rss <- function(x) var(x) * (length(x) - 1)

gsTree <- function(formula, input, model.control= list(minS = 20, minD = 5, error = c("deviance", "gini"))) {
  # Do not split if number of observations is less or equal to minS
  # Do not split if drop in rss is less or equal to minD percent
  # First we implement minS only;

    input <- as.data.frame(input)

    formula_terms <- terms(formula, data = input)
    # stop if formula has no response
    stopifnot(attr(formula_terms, 'response') == 1)

    Y_label <- as.character(formula_terms)[2]
    Y <- model.response(model.frame(formula, input))

    # X_label <- attr(formula_terms, "terms.labels")
    X <- as.data.frame(model.matrix(formula, input)[, -1])
    X_label <- colnames(X)

    if( is.factor(Y) ){
      stopifnot( ! is.factor(Y) ) # Classification trees aren't implemented yet.
    } else {
      minFUN <- rss
    }
    min_size <- model.control$minS

    split_along_predictor <- function(x, y) {
      if(is.factor(x)){
        splits <- levels(x)
        op <- '=='
      } else {
        splits <- sort(unique(x))
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
          rv <- c(error_left=minFUN(y[split_filter]), error_right=minFUN(y[!split_filter]))
        return(rv)
      })
      error_value <- sapply(split_list, sum)
      if(all(is.na(error_value))) return(NULL) # not possible to split without breaking min_size requirement
      i <- which.min(error_value)
      rv <- list(split = splits[i], op = op, error = split_list[[i]])
      return(rv)
    }

    bestNodeSplit <- function(node) {
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
          parent = i
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
        if(is.null(tmp)) {
          child$status <- 'L'  # no allowed splits
          child['best_split'] <- list(tmp)      # this is how one can assign NULL to a list element;
        }

        children[[j]] <- child
      }
      stopifnot((children[[1]]$obsN + children[[2]]$obsN) == node_to_split$obsN)
      return(children)
    }
    
    markNodeLeaf <- function(node){
      node$status <- 'L'
      node['best_split'] <- list(NULL)
      return(node)
    }
   
    # Tree to be grown
    rtree_ <- list()
    dtree_ <- NULL
    leaves_ <- NULL
    # error value of the model
    errorVal <- minFUN(Y)
    
    

    root_node <- list(
      nN = 1,                        # Node number
      status = 'S',
      parent = NULL,                 # Number of a parent Node
      split_conditions = "TRUE",       # List of split conditions to create this Node
      split_conditions_list = NULL,
      value = mean(Y),
      error = minFUN(Y),
      best_split = NULL,           # list with the following members:
      # split's feature and split's value, rss_l of left child node,
      # error_right of the right child node
      obsN = nrow(input)           # Number of observations for this Node
    )
    
    tmp <- bestNodeSplit(root_node)
    if(is.null(tmp)) { # no allowed splits
      stop("Tree can't be build under current conditions, try changing model.control argument")
    }
    root_node['best_split'] <- list(tmp)
    rtree_[[1]] <- root_node

    nodesN_to_split <- 1
    while(length(nodesN_to_split) > 0){

      cat("Nodes that can be split: ", paste(nodesN_to_split, collapse = ", "), "\n")

      best_split <- findBestNodeToSplit(nodesN_to_split)
      cat("Best split: ", paste(best_split, collapse = ", "), "\n")
      # update parent
      rtree_[[best_split]]$status<- 'P'
      nodeToSplit <- rtree_[[best_split]]
      deltaError <- with(nodeToSplit, error - sum(best_split$error))
      if(100 * deltaError / errorVal[length(errorVal)] < model.control$minD){
        for(i in nodesN_to_split) rtree_[[i]] <- markNodeLeaf(rtree_[[i]])
        break
      } else {
        children <- createSplit(best_split)
        rtree_ <- c(rtree_, children)
        errorVal <- c(errorVal, sum(sapply(rtree_, function(x){
          if(x$status == 'P') return(0) else return(x$error)
        })))
      }
      nodesN_to_split <- seq_len(length(rtree_))[sapply(rtree_, function(x) x$status == 'S')]

    }
    
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
       # nm <- paste(rtree_node$split_conditions_list[[length(rtree_node$split_conditions_list)]], collapse = " ")
        nm <- with(rtree_node, split_conditions[length(split_conditions)])
      } else {
        nm <- nodeDefnStr(rtree_node)
      }
      rv <- Node$new(nm)
      rv$parent_nN <- rtree_node$parent
      for(nm in setdiff(names(rtree_node), 'parent')) rv[[nm]] <- rtree_node[[nm]]
      rv
    }
    
    toDataTree <- function(tri){
      library(data.tree)
      
      flat_tree <- lapply(tri, toDataTreeNode)
      assembled_tree <- flat_tree[[1]]
      for(i in 2 : length(flat_tree)){
        flat_tree[[flat_tree[[i]]$parent_nN]]$AddChildNode(flat_tree[[i]])
      }
      return(assembled_tree)
    }
    
    dataTree <- function(){
      if( is.null(dtree_) ) dtree_ <<- toDataTree(rtree_)
      return(dtree_)
    }
    
    printTree <- function(){
      base::print(dataTree(), "obsN", "value", "error")
    }
    
    extractLeaves <- function(as.dataframe = FALSE){
      if(is.null(leaves_)){
        leaves_index <- seq_along(rtree_)[sapply(rtree_, function(x) x$status == 'L')]
        leaves <- lapply(leaves_index, function(x) {
          xn <- rtree_[[x]]
          xn[["defn"]] <- nodeDefnStr(xn)
          xn <- xn[c("nN", "split_conditions", "value", "error", "obsN", "defn")]
          xn
        })
        leaves_ <<- leaves
      }
      if(isTRUE(as.dataframe)){
        rv <- lapply(leaves_, function(x) {
          x <- within(x, split_conditions <- paste(split_conditions, collapse = " & "))
          x <- as.data.frame(x, stringsAsFactors = FALSE)
          row.names(x) <- NULL
          x
        })
        return(do.call(rbind, rv))
      }
      return(leaves_)
    }
        
    printLeaves <- function(){
      print(extractLeaves(as.dataframe = TRUE))
    }
    
    prediction <- function(X){
      leaves <- extractLeaves(as.dataframe = TRUE)
      Y <- vector("list", nrow(X))
      for(i in seq_len(nrow(leaves))){
        nodeF <- eval(str2lang(leaves[i, "split_conditions"]), envir=X)
        Y[nodeF] <- leaves[i, 'value']
      }
      return(unlist(Y))
    }
    
    plotTree <- function(){
      plot(dataTree())
    }
    
    
    rv <- list(
      tree = rtree_,
      errVal = errorVal,
      printTree = printTree,
      printLeaves = printLeaves,
      plotTree = plotTree,
      predict = prediction
    )
    return(rv)
  }

  #Example:
  library(MASS)
  library(tree)
  data(Boston)

  set.seed(1)
  train <- sample(seq_len(nrow(Boston)), nrow(Boston)/2)

  #foo <- gsTree(medv ~ ., input = Boston[train, ], model.control=list(minS = 5, minD = 0))
  foo <- gsTree(medv ~ ., input = Boston[train, ])
  Y <- foo$predict(Boston[-train, ])

  fit <- tree(medv ~ ., data = Boston[train, ])
  YY <- predict(fit, newdata = Boston[-train, ])
  
  data.frame(Y = Y, YY = YY, YYY = Boston[-train, 'medv'])