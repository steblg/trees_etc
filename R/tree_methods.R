
# Tree methods

information_gain <- function(xtree) {
  return(sapply(xtree, function(x) {
    rv <- with(x, {if (is.null(best_split)) 0 else sum(best_split$error)})
  }))
}

substitution_error <- function(xtree) {
  info_gain <- information_gain(xtree)
  N <- length(info_gain)
  sub_error <- vector(mode='numeric', length = N)
  for (i in N:1) {
    sub_error[i] <- with(xtree[[i]], {if (is.null(children_nN)) 0 else info_gain[i] - (sub_error[children_nN[1]] + sub_error[children_nN[2]])})
  }
  sub_error
}