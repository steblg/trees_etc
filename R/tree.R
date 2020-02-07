

# ToDo:
#   Handling of NA's. For now one has to deal with them (na.omit) prior to using the method;
#   Feature importance
#   Move to S3
#   Package namespace
#   Change "source" to a proper "package" way of loading files

source(file.path(getwd(), "R", "combinations.R"))
source(file.path(getwd(), "R", "basic_tree_defn.R"))
source(file.path(getwd(), "R", "basic_tree.R"))
source(file.path(getwd(), "R", "tree_methods.R"))
source(file.path(getwd(), "R", "complexity_pruning.R"))
source(file.path(getwd(), "R", "tree_pruning.R"))


#Example:
library(MASS)
# library(tree)
library(rpart)
library(rpart.plot)
data(Boston)

set.seed(1)
train <- sample(seq_len(nrow(Boston)), nrow(Boston)/2)

if (0) {
  # Compare simple tree building:
  reg_tree <- basic_tree(medv ~ ., input = Boston[train, ])
  Y <- predict_values(reg_tree, Boston[-train, ])
  # fit <- tree(medv ~ ., data = Boston[train, ]) # This tree corresponds to model.control = list(min_node_size = 5, min_sub_err = 10)
  fit <- rpart(medv ~ ., data = Boston[train, ], control = rpart.control(minbacket = 10, cp = .009)) # This should produce the same tree as "reg_tree"
  YY <- predict(fit, newdata = Boston[-train, ])

  YYY <- Boston[-train, 'medv']
  result <- data.frame(Y = Y, YY = YY, YYY = YYY)

  merr <- sum((Y - YYY)^2)/length(Y)
  terr <- sum((YY - YYY)^2)/length(YY)

  message(paste("My Err:", merr, "Prod Err:", terr, "\n"))

  stop()
}

if (1) {
  # Pruning
  
  # mfit <- basic_tree(medv ~ ., input = Boston[train, ], model.control = model_control(min_node_size = 10, min_sub_err = 1))
  # this mfit is approximately equal to the below rpart_fit
  rpart_fit <- rpart(medv ~ ., data = Boston[train, ], control = rpart.control(minbacket = 5, cp = .001))
  plot(rpart_fit)
  text(rpart_fit, pretty = FALSE, cex = .7)
  #cost-complexity pruning
  printcp(rpart_fit)
  opt <- which.min(rpart_fit$cptable[,"xerror"])
  #get its value
  cp <- rpart_fit$cptable[opt, "CP"]
  #prune tree
  rpart_fit_pruned <- prune(rpart_fit,cp)
  #plot tree
  plot(rpart_fit_pruned);
  text(rpart_fit_pruned, pretty = FALSE, cex = .7)
  
  prune_info <- cvtune(medv ~ ., input = Boston[train, ], k = 10, model.control = model_control(min_node_size = 5, min_sub_err = .5))
  print(as.data.frame(prune_info[c('betas', 'avg_error')]))
  ind <- 24
  print_tree(prune_info$tree, active = prune_info$prune_seq$info[[ind]]$active)
  prune_info_leaves <- sapply(prune_info$prune_seq$info, function(x) sum(x$leaves))
  prune_info_df <- cbind(as.data.frame(prune_info[c('betas', 'avg_error')]), complexity = prune_info_leaves)
  plot(avg_error ~ complexity, data = prune_info_df)
  
  mfit <- prune_info$tree
  print_tree(mfit)
  
  
  Y <- predict_values(prune_info$tree, Boston[-train, ], active = prune_info$prune_seq$info[[ind]]$active)
  YY <- predict(rpart_fit_pruned, Boston[-train, ])

  YYY <- Boston[-train, 'medv']
  result <- data.frame(Y = Y, YY = YY, YYY = YYY)

  merr <- sum((Y - YYY)^2)/length(Y)
  terr <- sum((YY - YYY)^2)/length(YY)

  message(paste("My Err:", merr, "Prod Err:", terr, "\n"))
  
  stop()
}
