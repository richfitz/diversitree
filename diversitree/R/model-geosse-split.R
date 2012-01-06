## GeoSSE model, by Emma Goldberg <eeg@uic.edu>

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

## 1: make
make.geosse.split <- function(tree, states, nodes, split.t,
                             sampling.f=NULL, strict=TRUE,
                             control=list()) {
  control <- check.control.ode(control)
  if ( control$backend == "CVODES" )
    stop("Cannot use CVODES backend with geosse.split")

  cache <- make.cache.geosse.split(tree, states, nodes, split.t, 
                                   sampling.f=sampling.f,
                                   strict=strict)
  cache$control <- check.control.split(control)

  branches.main <- make.branches.geosse(cache, control)
  branches.aux <- make.branches.aux.geosse(cache, control)
  branches <- make.branches.split(cache, branches.main, branches.aux)
  initial.conditions <-
    make.initial.conditions.split(cache, initial.conditions.geosse)

  n.part <- cache$n.part
  
  ## (note: condition.surv=TRUE in bisse)
  ll <- function(pars, condition.surv=FALSE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.par.multipart(pars, n.part, 7)
    pars.n <- unlist(pars)
    if ( any(pars.n < 0) || any(!is.finite(pars.n)) )
      stop("All parameters must be nonnegative")

    if ( !is.null(root.p) && root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    pars.root <- pars[[cache$group.nodes[cache$root]]]
    if ( cache$control$caching.branches )
      caching.branches.set.pars(pars, branches)

    ans <- all.branches.matrix(pars, cache, initial.conditions,
                               branches)

    vals <- ans$init[,cache$root]
    root.p <- root.p.geosse(vals, pars.root, root, root.p)
    loglik <- root.geosse(vals, pars.root, ans$lq, condition.surv, root.p)

    if (intermediates) {
        ans$root.p <- root.p
        attr(loglik, "intermediates") <- ans
        attr(loglik, "vals") <- vals
    }

    loglik
  }
 
  class(ll) <- c("geosse.split", "geosse", "function")
  attr(ll, "n.part") <- cache$n.part
  ll
}

## 2: print
print.geosse.split <- function(x, ...) {
  cat("GeoSSE(split) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.geosse.split <- function(x, n.part=attr(x, "n.part"), ...) {
  argnames.twopart(x, argnames.geosse(NULL), n.part)
}
`argnames<-.geosse.split` <- function(x, value) {
  n.part <- attr(x, "n.part")
  argnames.twopart.set(x, value, 7, n.part)
}

## 4: find.mle
find.mle.geosse.split <- function(func, x.init, method, fail.value=NA,
                                 ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method,
             class.append="fit.mle.geosse.split")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.geosse.split <- function(tree, states, nodes, split.t,
                                    sampling.f, strict) {
  cache <- make.cache.geosse(tree, states, NULL, sampling.f, NULL,
                             strict)
  cache <- make.cache.split(tree, cache, nodes, split.t)
  cache$sampling.f <- check.sampling.f.split(sampling.f, 3, cache$n.part)
  cache$aux.i <- 1:3
  cache$y <- make.initial.tip.xxsse.split(cache)  
  cache
}

## 7: initial.conditions: from geosse

## 8: branches: from geosse.  However the 'branches.aux' function is
## required to compute the E0, E1, E2 values after a partition.
make.branches.aux.geosse <- function(cache, control) {
  idx.e <- 1:3
  y <- lapply(cache$sampling.f, function(x) c(1-x, rep(1, 3)))
  n <- length(y)
  branches <- make.branches.geosse(cache, control)

  function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[[2]][idx.e,,drop=FALSE]
  }
}
