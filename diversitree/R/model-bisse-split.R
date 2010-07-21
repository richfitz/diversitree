## Models should provide:
##   1. make
##   2. print
##   3. argnames / argnames<-
##   4. find.mle
## Generally, make will require:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches
##   9. branches.unresolved

## 1: make
make.bisse.split <- function(tree, states, nodes, split.t,
                             unresolved=NULL, sampling.f=NULL,
                             nt.extra=10, safe=FALSE) {
  cache <- make.cache.bisse.split(tree, states, nodes, split.t,
                                  unresolved, sampling.f, nt.extra)
  branches <- make.branches.bisse(safe)
  branches.aux <- make.branches.aux.bisse(cache$sampling.f, safe)
  ll <- function(pars, ...)
    ll.bisse.split(cache, pars, branches, branches.aux, ...)
  class(ll) <- c("bisse.split", "function")
  attr(ll, "n.part") <- cache$n.part
  ll
}

## 2: print
print.bisse.split <- function(x, ...) {
  cat("BiSSE(split) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.bisse.split <- function(x, ...) {
  obj <- attr(x, "argnames")
  n <- attr(x, "n")
  if ( is.null(obj) )
    obj <- list(base=c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10"),
                levels=seq_len(n))

  paste(obj$base, rep(obj$levels, each=6), sep=".")
}
`argnames<-.bisse.split` <- function(x, value) {
  n <- attr(x, "n")
  if ( !is.list(value) || length(value) != 2 )
    stop("'value' must be a list of length 2")
  if ( length(value[[1]]) != 6 || length(value[[2]]) != n )
    stop(sprintf("value's elements must be of length 6, %d", n))

  names(value) <- c("base", "levels")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.bisse.split <- function(func, x.init, method, fail.value=NA,
                                 ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method,
             class.append="fit.mle.bisse.split")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.bisse.split <- function(tree, states, nodes, split.t,
                                   unresolved=NULL, sampling.f=NULL,
                                   nt.extra=10) {
  ## 1: tree
  tree <- check.tree(tree)

  ## 2: states:
  states <- check.states(tree, states)

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  ## TODO: a sampling.f of length 2 can be expanded manually.
  n <- length(nodes) + 1 # +1 for base group
  if ( is.list(sampling.f) ) {
    check.sampling.f(unlist(sampling.f), 2*n)
  } else {
    sampling.f <- check.sampling.f(sampling.f, 2 * n)
    sampling.f <- matrix.to.list(matrix(sampling.f, n, 2, TRUE))
  }

  cache <- make.cache.split(tree, nodes, split.t)

  for ( i in seq_along(cache$cache) ) {
    x <- cache$cache[[i]]
    x$tip.state  <- states[x$tip.label]
    x$sampling.f <- sampling.f[[i]]

    if ( !is.null(unresolved) ) {
      ok <- unresolved$tip.label %in% x$tip.label
      if ( any(ok) ) {
        x$unresolved <-
          check.unresolved(x, unresolved[ok,], nt.extra)
        x$tips <- x$tips[-x$unresolved$i]
        x$tip.label <- x$tip.label[-x$unresolved$i]
        x$tip.state <- x$tip.state[-x$unresolved$i]
      }
    }
    
    x$y <- initial.tip.bisse(x)
    cache$cache[[i]] <- x
  }

  cache$sampling.f <- sampling.f
  cache$aux.i <- 1:2
  cache
}

ll.bisse.split <- function(cache, pars, branches, branches.aux,
                           condition.surv=TRUE, root=ROOT.OBS,
                           root.p=NULL, intermediates=FALSE) {
  n.part <- cache$n.part
  if ( is.matrix(pars) ) {
    if ( nrow(pars) != n.part )
      stop(sprintf("Expected %d parameter sets", n.part))
    if ( ncol(pars) != 6 )
      stop("Expected 6 parameters in each set")
    pars <- matrix.to.list(pars)
  } else if ( is.list(pars) ) {
    if ( length(pars) != n.part )
      stop(sprintf("Expected %d parameter sets", n.part))
    if ( !all(unlist(lapply(pars, length)) == 6) )
      stop("Expected 6 parameters in each set")
  } else {
    if ( length(pars) != n.part * 6 )
      stop(sprintf("Expected %d parameters", n.part * 6))
    pars <- matrix.to.list(matrix(pars, n.part, 6, TRUE))
  }
  pars.n <- unlist(pars)
  if ( any(pars.n < 0) || any(!is.finite(pars.n)) )
    return(-Inf)

  for ( i in seq_len(n.part) ) {
    unresolved.i <- cache$cache[[i]]$unresolved
    if ( !is.null(unresolved.i) )
      cache$cache[[i]]$preset <-
        branches.unresolved.bisse(pars[[i]], unresolved.i)
  }

  ans <- all.branches.split(pars, cache, initial.conditions.bisse,
                            branches, branches.aux)

  vars <- ans[[1]]$base
  lq <- unlist(lapply(ans, "[[", "lq"))
  pars.root <- pars[[1]]

  root.p <- root.p.xxsse(vars, pars.root, root, root.p)
  loglik <- root.xxsse(vars, pars.root, lq, condition.surv, root.p)

  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
}

## 7: initial.conditions: from bisse

## 8: branches: from bisse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.

## TODO: this would be nicer if it did not compute the Ds at all.
## Also, it can just use the non-clever function as I do not need the
## log compensation worked out.
make.branches.aux.bisse <- function(sampling.f, safe=FALSE) {
  y <- lapply(sampling.f, function(x) c(1-x, 1, 1))
  branches <- make.branches.bisse(safe)
  n <- length(y)
  branches.aux.bisse <- function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[,2:3,drop=FALSE]
  }
}


## 9: branches.unresolved: from bisse
