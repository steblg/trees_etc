
combinations_old <- function(x, k) {
	# Unique combinations of length 'k', up to permutations, 
	# from elements of 'x'.
	# If 'x' is a list, combinations are returned
	# as vectors of indices to the list.
	
	if ( is.list(x) ) x <- seq_along(x)
	stopifnot(is.vector(x) && is.atomic(x))
	stopifnot(!any(duplicated(x)))

	n <- length(x)
	stopifnot(k <= n)
	
	if (k == 1) return(as.list(x))
	if (n == k) return(list(x))
	
	rv <- c(lapply(combinations_old(x[-1], k - 1), function(y) c(x[1], y)), combinations_old(x[-1], k))
	rv
}

combinations <- function(x, k) {
	# Unique combinations of length 'k', up to permutations, 
	# from elements of 'x'.
  # Elements of 'x' must be unique
	# If 'x' is a list, combinations are returned
	# as vectors of indices to the list.
	
	if ( is.list(x) ) x <- seq_along(x)
	stopifnot(is.vector(x) && is.atomic(x))
	stopifnot(!any(duplicated(x)))
	
	n <- length(x)
	stopifnot(k <= n)
	
	if (k == 1) return(as.list(x))
	if (n == k) return(list(x))
	
	rv <- lapply(1:(n - k + 1), function(i) {
		xx <- x[(i + 1) : n]
		lapply(combinations(xx, k - 1), function(y) c(x[i], y))
	})
	rv <- do.call(c, rv)
	rv
}

sample_space <- function(x, complimentary = FALSE) {
  # All unique combinations of all possible lengths, up
  # to permutations, from elements of 'x'
  # Elements of 'x' must be unique. 
  # Since length of 'x' is finite and all 'x' elements are unique,
  # every unique combination "C_1" of elements of 'x' has unique
  # combination "C_2" that is complimentary, i.e. C_2 == !C_1.
  # By default, these complimentary combinations are not generated.
  # To include them, set "complimentary" argument to TRUE
	stopifnot(is.vector(x) && is.atomic(x))
	stopifnot(!any(duplicated(x)))
	
	n <- length(x)
	if (isTRUE(complimentary)) maxk <- n else maxk <- n %/% 2
	
	rv <- lapply(1:maxk, function(i) combinations(x, i))
	rv <- do.call(c, rv)
	rv
}
