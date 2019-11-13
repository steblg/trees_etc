

# ToDo:
#   1. Handling of NA's. For now one has to deal with them (na.omit) prior to using the method;
#   2. Feature importance
#   3. Complexity pruning w cross-validation

#   5. Move to S3
#   6. Package namespace
#   7. Change "source" to a proper "package" way of loading files

source(file.path(getwd(), "R", "basic_tree.R"))
source(file.path(getwd(), "R", "tree_methods.R"))
source(file.path(getwd(), "R", "complexity_pruning.R"))
source(file.path(getwd(), "R", "tree_pruning.R"))


#Example:
library(MASS)
library(tree)
data(Boston)

set.seed(1)
train <- sample(seq_len(nrow(Boston)), nrow(Boston)/2)

reg_tree <- basic_tree(medv ~ ., input = Boston[train, ], model.control=list(minS = 5, minD = 5))
# reg_tree <- basic_tree(medv ~ ., input = Boston[train, ])
Y <- predict_values(reg_tree, Boston[-train, ])

fit <- tree(medv ~ ., data = Boston[train, ])
YY <- predict(fit, newdata = Boston[-train, ])
YYY <- Boston[-train, 'medv']
result <- data.frame(Y = Y, YY = YY, YYY = YYY)

merr <- sum((Y - YYY)^2)/length(Y)
terr <- sum((YY - YYY)^2)/length(YY)

prune_info <- cvtune(medv ~ ., input = Boston[train, ], model.control=list(minS = 5, minD = 5))
print_tree(prune_info$tree, active = prune_info$prune_seq$active[[7]])

print_tree(reg_tree, cc_prune_tree(reg_tree, alpha = 200)[[1]])