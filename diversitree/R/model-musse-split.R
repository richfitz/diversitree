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
make.musse.split <- function(tree, states, k, nodes, split.t,
                             sampling.f=NULL, strict=TRUE,
                             control=list()) {
  control <- check.control.ode(control)
  if ( control$backend == "CVODES" )
    stop("Cannot use CVODES backend with musse.split")

  cache <- make.cache.musse.split(tree, states, k, nodes, split.t,
                                  sampling.f, strict)

  branches <- make.branches.musse(cache, control)
  branches.aux <- make.branches.aux.musse(cache, control)

  k <- cache$k
  n.part <- cache$n.part
  np <- k * (k + 1)

  f.pars <- make.musse.pars(k)  
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.par.multipart(pars, n.part, np)
    pars.n <- unlist(pars)
    if ( any(pars.n < 0) || any(!is.finite(pars.n)) )
      return(-Inf)

    ## Set up all the Q matrices:
    for ( i in seq_len(n.part) )
      pars[[i]] <- f.pars(pars[[i]])

    ans <- all.branches.split(pars, cache, initial.conditions.musse,
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
 
  class(ll) <- c("musse.split", "musse", "function")
  attr(ll, "n.part") <- cache$n.part
  attr(ll, "k")      <- k
  ll
}

## 2: print
print.musse.split <- function(x, ...) {
  cat("Musse(split) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.musse.split <- function(x, k=attr(x, "k"),
                                 n.part=attr(x, "n.part"), ...) {
  argnames.twopart(x, argnames.musse(NULL, k), n.part)
}
`argnames<-.musse.split` <- function(x, value) {
  n.part <- attr(x, "n.part")
  k <- attr(x, "k")
  argnames.twopart.set(x, value, k * (k + 1), n.part)
}

## 4: find.mle
find.mle.musse.split <- function(func, x.init, method, fail.value=NA,
                                 ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method,
             class.append="fit.mle.musse.split")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.musse.split <- function(tree, states, k, nodes, split.t,
                                   sampling.f, strict) {
  ## 1: tree
  tree <- check.tree(tree, node.labels=TRUE)

  ## 2: states:
  states <- check.states(tree, states,
                         strict=strict, strict.vals=1:k)

  n <- length(nodes) + 1 # +1 for base group
  sampling.f <- check.sampling.f.split(sampling.f, k, n)

  cache <- make.cache.split(tree, nodes, split.t)

  for ( i in seq_along(cache$cache) ) {
    x <- cache$cache[[i]]
    x$ny <- 2*k
    x$k <- k
    x$tip.state  <- states[x$tip.label]
    x$sampling.f <- sampling.f[[i]]
    x$y <- initial.tip.musse(x)
    cache$cache[[i]] <- x
  }

  cache$sampling.f <- sampling.f
  cache$aux.i <- 1:k
  cache$k <- k

  cache
}

## 7: initial.conditions: from musse

## 8: branches: from musse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.
make.branches.aux.musse <- function(cache, control) {
  k <- cache$k
  idx.e <- seq_len(k)
  y <- lapply(cache$sampling.f, function(x) c(1-x, rep(1, k)))
  n <- length(y)
  branches <- make.branches.musse(cache, control)

  function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[[2]][idx.e,,drop=FALSE]
  }
}
