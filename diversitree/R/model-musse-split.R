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
  backend <- control$backend

  cache <- make.cache.musse.split(tree, states, k, nodes, split.t,
                                  sampling.f, strict)
  cache$control <- check.control.split(control)
  n.part <- cache$n.part
 
  if ( backend == "CVODES" )
    stop("This will take some work...")
  else {
    branches.main <- make.branches.musse(cache, control)
    branches.aux <- make.branches.musse.aux(cache, control)
    branches <- make.branches.split(cache, branches.main, branches.aux)    
    initial.conditions <-
      make.initial.conditions.split(cache, initial.conditions.musse)
  }

  n.part <- cache$n.part
  k <- cache$k
  np <- k * (k + 1)

  f.pars <- make.musse.pars(k)  
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.par.musse.split(pars, n.part, k)
    for ( i in seq_len(n.part) )
      pars[[i]] <- f.pars(pars[[i]])

    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    if ( backend == "CVODES" )
      stop("Not done yet...")
    else 
      ll.xxsse.split(pars, cache, initial.conditions, branches,
                     condition.surv, root, root.p, intermediates)
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
  cache <- make.cache.musse(tree, states, k, NULL, strict)
  cache <- make.cache.split(tree, cache, nodes, split.t)
  cache$sampling.f <- check.sampling.f.split(sampling.f, k, cache$n.part)
  cache$aux.i <- seq_len(k)
  cache$y <- make.initial.tip.xxsse.split(cache)  
  cache
}

## 7: initial.conditions: from musse

## 8: branches: from musse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.
make.branches.musse.aux <- function(cache, control) {
  k <- cache$k
  np <- as.integer(k * (k + 2))
  neq <- as.integer(k)
  comp.idx <- integer(0)
  branches <- make.ode.branches("musse_aux", "diversitree", neq, np,
                                comp.idx, control)
  y <- lapply(cache$sampling.f, function(x) 1-x)
  n <- length(y)  

  function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[[2]]
  }
}
