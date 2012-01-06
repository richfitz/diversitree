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
                             nt.extra=10, strict=TRUE,
                             control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend
  if ( backend == "CVODES" && !is.null(cache$unresolved) )
    stop("Cannot yet use CVODES backend with unresolved clades")

  cache <- make.cache.bisse.split(tree, states, nodes, split.t,
                                  unresolved, sampling.f, nt.extra,
                                  strict)
  cache$control <- check.control.split(control)
  
  if ( backend == "CVODES" )
    stop("This will take some work...")
  else {
    branches.main <- make.branches.bisse(cache, control)
    branches.aux <- make.branches.bisse.aux(cache, control)
    branches <- make.branches.split(cache, branches.main, branches.aux)
    initial.conditions <-
      make.initial.conditions.split(cache, initial.conditions.bisse)
  }

  n.part <- cache$n.part
  unresolved <- cache$unresolved
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.par.bisse.split(pars, n.part)

    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    if ( !is.null(unresolved) )
      cache$preset <-
        branches.unresolved.bisse.split(pars, unresolved)

    if ( backend == "CVODES" )
      stop("Not yet done...")
    else
      ll.xxsse.split(pars, cache, initial.conditions, branches,
                     condition.surv, root, root.p, intermediates)
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
                                   nt.extra=10, strict=TRUE) {
  cache <- make.cache.bisse(tree, states, unresolved, NULL, nt.extra,
                            strict)
  cache <- make.cache.split(tree, cache, nodes, split.t)
  cache$sampling.f <- check.sampling.f.split(sampling.f, 2, cache$n.part)
  cache$aux.i <- 1:2

  unresolved <- cache$unresolved
  n.part <- cache$n.part

  if ( !is.null(unresolved) ) {
    ## This ensures that the calculations should be slightly more
    ## identical by running out to the same number of species as for
    ## the non-split version.  Where one group mostly has smaller
    ## clades, this will slow the calculations down more than needed
    ## (will drop this after things settle down).
    Nc.tot <- max(unresolved$Nc) + unresolved$nt.extra

    grp <- cache$group.branches[unresolved$target]
    tmp <- data.frame(unresolved[names(unresolved) != "nt.extra"],
                                   stringsAsFactors=FALSE)
    unresolved.split <- vector("list", n.part)
    
    for ( i in seq_len(n.part)[tabulate(grp, n.part) > 0] ) {
      j <- grp == i
      unresolved.split[[i]] <- as.list(tmp[j,])
      unresolved.split[[i]]$nt.extra <-
        Nc.tot - max(unresolved.split[[i]]$Nc)
    }

    cache$unresolved <- unresolved.split
  }

  cache$y <- make.initial.tip.xxsse.split(cache)
  
  cache
}

## 7: initial.conditions: from bisse

## 8: branches: from bisse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.
make.branches.bisse.aux <- function(cache, control) {
  neq <- 2L
  np <- 6L
  comp.idx <- integer(0)
  branches <- make.ode.branches("bisse_aux", "diversitree", neq, np,
                                comp.idx, control)
  y <- lapply(cache$sampling.f, function(x) 1-x)
  n <- length(y)  

  function(i, len, pars) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[[2]]
  }
}

branches.unresolved.bisse.split <- function(pars, unresolved) {
  ans <- list(target=integer(0), lq=numeric(0),
                 base=matrix(NA, 4, 0))
  for ( i in seq_along(pars) ) {
    unresolved.i <- unresolved[[i]]
    if ( !is.null(unresolved.i) ) {
      tmp <- branches.unresolved.bisse(pars[[i]], unresolved.i)
      ans$target <- c(ans$target, tmp$target)
      ans$lq <- c(ans$lq, tmp$lq)
      ans$base <- cbind(ans$base, tmp$base)
    }
  }
  ans
}

check.par.musse.split <- function(pars, n.part, k) {
  pars <- check.par.multipart(pars, n.part, k*(k+1))
  pars.n <- unlist(pars) # for checking only...
  if ( any(pars.n < 0) || any(!is.finite(pars.n)) )
    stop("All parameters must be nonnegative")
  pars
}

check.par.bisse.split <- function(pars, n.part)
  check.par.musse.split(pars, n.part, 2)
