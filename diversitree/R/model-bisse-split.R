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
  ll <- function(pars, ...)
    ll.bisse.split(cache, pars, branches, ...)
  class(ll) <- c("bisse.split", "function")
  attr(ll, "n") <- cache$n
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

  matrix(paste(rep(obj$base, each=n), obj$levels, sep="."), n, 6)
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
                                   nt.extra=10, safe=FALSE) {
  ## 1: tree
  tree <- check.tree(tree)

  ## 2: states:
  states <- check.states(tree, states)

  ## 3: A large amount of unstreamlined checking of sampling.f and
  ## unresolved - try to abstract this to some degree.

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  n <- length(nodes) + 1 # +1 for base group
  sampling.f <- check.sampling.f(sampling.f, 2 * n)
  sampling.f <- matrix(sampling.f, n, 2)

  ## This has the poor effect of not working correctly to create
  ## single branch partitions.  I should be careful about that.
  subtrees <- split.phylo(phy, nodes, split.t=split.t)

  ## Next, process the unresolved clade information for the different
  ## trees, and build the

  if ( !is.null(unresolved) ) {
    ret <- vector("list", n)
    for ( i in seq_len(n) ) {
      x <- subtrees[[i]]
      ok <- (unresolved$tip.label %in%
             setdiff(x$tip.label, names(x$daughters)))
      if ( any(ok) )
        ret[[i]] <- check.unresolved(x, unresolved[ok,], nt.extra)
    }
    unresolved <- ret
  }
  
  cache <- list()
  cache$cache <- vector("list", n)
  for ( i in seq_len(n) )
    cache$cache[[i]] <-
      mcbs.one(subtrees[[i]], states, unresolved[[i]], sampling.f[i,])
  
  cache$daughters <- lapply(subtrees, "[[", "daughters")
  cache$parents <- sapply(subtrees, "[[", "parent")
  cache$order <- dt.split.order(cache$daughters, cache$parents)
  cache$n <- length(subtrees)

  cache
}

mcbs.one <- function(tree.sub, states, unresolved, sampling.f) {
  ## This is essentially identical to the code in make.cache.bisse()
  cache <- make.cache(tree.sub)
  cache$tip.state  <- states[tree.sub$tip.label]
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f
  cache$y <- initial.tip.bisse(cache)

  ## This adds the daughter information:
  cache$daughters <- tree.sub$daughters
  n <- length(tree.sub$daughters)
  cache$daughters.i <- seq_len(n) + 3
  if ( n > 0 ) {
    cache$y$y <- rbind(cache$y$y, matrix(NA, n, 4))
    cache$y$i[names(tree.sub$daughters)] <- cache$daughters.i
    cache$y$types <- c(cache$y$types, cache$daughters.i)
  }

  cache$trailing <- tree.sub$trailing
  cache$trailing.t0 <- max(cache$depth)
  
  cache
}


ll.bisse.split <- function(cache, pars, branches, condition.surv=TRUE,
                           root=ROOT.OBS, root.p=NULL,
                           intermediates=FALSE) {
  if ( intermediates )
    .NotYetUsed("intermediates")

  n <- cache$n
  
  if ( !is.matrix(pars) || all(dim(pars) != c(n, 6)) )
    stop(sprintf("Expected a %d x 6 matrix of parameters", n))
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)

  res <- matrix(NA, n, 5)
  obj <- vector("list", n)
  for ( i in cache$order ) {
    x <- cache$cache[[i]]

    if ( length(x$daughters) > 0 ) {
      ## Rather than take the E values from the daughter lineage, I
      ## need to recompute them to use the new parameters at this
      ## point.  This means that I have to run a new branch down with
      ## the parent set and parent sampling.f (I am using the "no
      ## data" initial condition here).  The D initial conditions come
      ## from the previous case.
      e0 <- branches(x$y$y[3,], x$depth[names(x$daughters)],
                     pars[i,], 0)[,2:3,drop=FALSE]
      d0 <- res[x$daughters,4:5,drop=FALSE]
      x$y$y[x$daughters.i,] <- cbind(e0, d0)
    }

    obj[[i]] <- all.branches(pars[i,], x, initial.conditions.bisse,
                             branches, branches.unresolved.bisse)

    ## Add trailing edge
    if ( !is.na(cache$parent[i]) )
      res[i,] <- branches(obj[[i]]$init[x$root,], x$trailing,
                          pars[i,], x$trailing.t0)
    else 
      res[i,] <- c(0, obj[[i]]$init[x$root,])

    res[i,1] <- res[i,1] + sum(obj[[i]]$lq)
  }
  
  root.p <- root.p.xxsse(res[i,-1], pars[i,], root, root.p)
  root.xxsse(res[i,-1], pars[i,], res[,1], condition.surv, root.p)
}

## 7: initial.conditions: from bisse
## 8: branches: from bisse
## 9: branches.unresolved: from bisse
