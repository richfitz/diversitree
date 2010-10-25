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
                             sampling.f=NULL, strict=TRUE, safe=FALSE) {
  cache <- make.cache.musse.split(tree, states, k, nodes, split.t,
                                  sampling.f, strict)
  branches <- make.branches.musse(k, safe)
  branches.aux <- make.branches.aux.musse(k, cache$sampling.f, safe)
  ll <- function(pars, ...)
    ll.musse.split(cache, pars, branches, branches.aux, ...)
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
    x$tip.state  <- states[x$tip.label]
    x$sampling.f <- sampling.f[[i]]
    x$k <- k
    x$y <- initial.tip.musse(x)
    cache$cache[[i]] <- x
  }

  cache$sampling.f <- sampling.f
  cache$aux.i <- 1:k
  cache$k <- k

  cache
}

ll.musse.split <- function(cache, pars, branches, branches.aux,
                           condition.surv=TRUE, root=ROOT.OBS,
                           root.p=NULL, intermediates=FALSE) {
  n.part <- cache$n.part
  k <- cache$k

  pars <- check.par.multipart(pars, n.part, k * (k + 1))

  pars.n <- unlist(pars)
  if ( any(pars.n < 0) || any(!is.finite(pars.n)) )
    return(-Inf)

  ans <- all.branches.split(pars, cache, initial.conditions.musse,
                            branches, branches.aux)

  vals <- ans[[1]]$base
  lq <- unlist(lapply(ans, "[[", "lq"))

  pars.root <- pars[[1]]
  root.p <- root.p.xxsse(vals, pars.root, root, root.p)
  loglik <- root.xxsse(vals, pars.root, lq, condition.surv, root.p)

  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
}

## 7: initial.conditions: from musse

## 8: branches: from musse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.

## TODO: this would be nicer if it did not compute the Ds at all.
## Also, it can just use the non-clever function as I do not need the
## log compensation worked out.
make.branches.aux.musse <- function(k, sampling.f, safe=FALSE) {
  y <- lapply(sampling.f, function(x) c(1-x, rep(1, k)))
  branches <- make.branches.musse(k, safe)
  n <- length(y)
  j <- 2:(k+1)
  branches.aux.musse <- function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[,j,drop=FALSE]
  }
}

