## Constrain an argument to another argument
constrain.par <- function(f, rel) {
  f <- match.fun(f)
  rel.i <- resolve.constraint(rel)

  if ( inherits(f, c("fixed", "constrained")) )
    stop("Cannot (yet) constrain a constrained function")
  else
    free <- is.na(rel)

  g <- function(x, ...)
    f(x[rel.i], ...)
  class(g) <- c(class(f), "constrained")
  g
}

## Fix an argument to a particular value.
fix.par <- function(f, rel) {
  f <- match.fun(f)

  if ( inherits(f, c("fixed", "constrained")) )
    ## free <- environment(f)$free & is.na(rel)
    stop("Cannot (yet) constrain a constrained function")    
  else
    free <- is.na(rel)

  g <- function(x, ...)
    f(spread(x, rel), ...)
  class(g) <- c(class(f), "fixed")
  g
}

spread <- function(x, rel) {
  idx <- is.na(rel)
  j <- which(!idx)
  y <- x[match(seq_along(rel), which(idx))]
  y[j] <- rel[j]
  y
}

resolve.constraint <- function(rel) {
  i <- is.na(rel)
  rel[i] <- seq_len(sum(i))
  rel[!i] <- rel[rel[!i]]
  rel
}

## Sensible starting point for optimisation
starting.point <- function(tree, q.div=5) {
 fit <- suppressWarnings(birthdeath(tree))
 r <- fit$para[2]
 e <- fit$para[1]
 p <- rep(c(r/(1-e), r*e/(1-e), r/q.div), each=2)
 names(p) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
 p
}

## Protect a function from failure.  This is a little like
## tryQuietly, but we will also ensure finiteness.
protect <- function(f, fail.value, finite=TRUE) {
  function(...) {
    ret <- trySilent(f(...))
    if ( inherits(ret, "try-error") ||
         (finite && (is.na(ret) || !is.finite(ret))) )
      fail.value
    else
      ret
  }
}

big.brother <- function(f, interval=1) {
  .x.eval <- list()
  .y.eval <- list()
  function(x, ...) {
    i <- length(.x.eval) + 1
    if ( i %% interval == 0 )
      cat(sprintf("[%s]", paste(formatC(x, 5), collapse=", ")))
    else
      cat(".")
    .x.eval[[i]] <<- x
    .y.eval[[i]] <<- ans <- f(x, ...)
    if ( i %% interval == 0 )
      cat(sprintf("\t -> %2.5f\n", ans))
    ans
  }
}
