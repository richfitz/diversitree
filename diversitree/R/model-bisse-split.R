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
                             nt.extra=10, control=list()) {
  control <- check.control.ode(control)
  if ( control$backend == "CVODES" )
    stop("Cannot use CVODES backend with bisse.split")

  cache <- make.cache.bisse.split(tree, states, nodes, split.t,
                                  unresolved, sampling.f, nt.extra)

  branches <- make.branches.bisse(cache, control)
  branches.aux <- make.branches.aux.bisse(cache, control)

  n.part <- cache$n.part

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.par.multipart(pars, n.part, 6)
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
                              branches, branches.aux, FALSE)

    vals <- ans[[1]]$base
    lq <- unlist(lapply(ans, "[[", "lq"))

    pars.root <- pars[[1]]
    root.p <- root.p.xxsse(vals, pars.root, root, root.p)
    loglik <- root.xxsse(vals, pars.root, lq, condition.surv, root.p)

    if ( intermediates ) {
      ans$root.p <- root.p
      attr(loglik, "intermediates") <- ans
      attr(loglik, "vals") <- vals
    }
    loglik
  }

  class(ll) <- c("bisse.split", "bisse", "function")
  attr(ll, "n.part") <- cache$n.part
  ll
}

## 2: print
print.bisse.split <- function(x, ...) {
  cat("BiSSE(split) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.bisse.split <- function(x, n.part=attr(x, "n.part"), ...) {
  argnames.twopart(x, argnames.bisse(NULL), n.part)
}
`argnames<-.bisse.split` <- function(x, value) {
  n.part <- attr(x, "n.part")
  argnames.twopart.set(x, value, 6, n.part)
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
  tree <- check.tree(tree, node.labels=TRUE)

  ## 2: states:
  states <- check.states(tree, states)

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  n <- length(nodes) + 1 # +1 for base group
  sampling.f <- check.sampling.f.split(sampling.f, 2, n)

  cache <- make.cache.split(tree, nodes, split.t)

  for ( i in seq_along(cache$cache) ) {
    x <- cache$cache[[i]]
    x$ny <- 4L
    x$k <- 2L
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

    if ( x$n.tip > 0 )
      x$y <- initial.tip.bisse(x)
    cache$cache[[i]] <- x
  }

  cache$sampling.f <- sampling.f
  cache$aux.i <- 1:2
  cache
}

## 7: initial.conditions: from bisse

## 8: branches: from bisse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.
make.branches.aux.bisse <- function(cache, control) {
  idx.e <- 1:2
  y <- lapply(cache$sampling.f, function(x) c(1-x, 1, 1))
  n <- length(y)
  branches <- make.branches.bisse(cache, control)

  function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[[2]][idx.e,,drop=FALSE]
  }
}
