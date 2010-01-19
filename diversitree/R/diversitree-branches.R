## This is the general calculation code for diversitree.  Many models
## require that calculations occur down a branch, and the results of
## these calculations are used as initial conditions for the next
## branch.  There are several different implemented models using this
## approach:
##   - BiSSE (model-bisse.R)
##   - MuSSE (multi-state BiSSE - model-musse.R)
##   - BiSSEtd (time-dependent BiSSE - model-bissetd.R)
##   - QuaSSE (MOL integrator only)
##   - BD2 (HIV model)

## I have abstracted most of the calculations in a (hopefully) general
## way.  The function 'all.branches' does most of the calculations;
## this looks after all the book-keeping as calculations proceed down
## a tree.  Models will generally use this as their main function, but
## with additional root state calculations after this calculation has
## finished (see model-bisse.R for the canonical example of this).

## all.branches takes arguments
##   pars, cache, initial.conditions, branches, branches.unresolved
## These are
##
##   pars: parameters for the model.  No checks are done at all on the
##   size or contents of this (all.branches knows nothing about what
##   is appropriate)
##
##   cache: A cache object produced by make.cache(), possibly
##   augmented with additional information.
##
##   initial.conditions: This calculates the initial conditions for an
##   *internal branch*, given the values of variables at the base of
##   its two daughters (this is not the initial conditions for the
##   tips, which are computed elsewhere - see make.cache).  This
##   function takes the arguments:
##     init, pars, is.root
##   where
##     'init': a two-row, npar column matrix of initial conditions
##     corresponding to the base conditions of the two daugher
##     branches.
##     'pars': the parameters as given into all.branches
##     'is.root': boolean indicating if the node is the root or not
##   This function must return a vector of length npar+1.  The first
##   of these is the log-compensation value (zero if none applied),
##   and the rest are the new initial conditions.  all.branches checks
##   to make sure that the log-compensation value is finite as a way
##   of making sure that the calculations succeeded.
##   initial.conditions may produce informative errors instead, or
##   warnings to inform of this failure.
##
##   branches: This function calculates values of variables at the
##   base of a branch, given that branchs' initial conditions.  It
##   takes arguments:
##     y, len, pars, t0
##   where
##     'y': a vector of initial conditions (length npar)
##     'len': the length of the branch
##     'pars': the parameters as passed to all.branches
##     't0': the time at the tip of the branch
##   Note that the base of the branch is at time t0 + len.  This
##   function must return a vector of length 'npar', being the
##   variables at the base of the branch.
##
##   'branches' may be able to deal with multiple branches at the
##   present (this is the case in BiSSE, where because there are only
##   three possible initial conditions, variables required for short
##   branches are often already computed on longer branches.  If
##   presented with a vector of 'len' at 'depth=0', branches must
##   return a matrix of variables of dimensions length(len) x npar.
##   Whether this is required or not depends on the initial conditions
##   produced by make.cache - see the documentation there and also for
##   the bisse initial.condition and branches functions.
##
##   branches.unresolved: This function computes the values of
##   variables for terminal unresolved clades.  It takes arguments
##     pars, len, unresolved
##   where
##     'pars': the parameters as passed to all.branches
##     'len': the ...
##   TODO: documentation coming - see source!
all.branches <- function(pars, cache, initial.conditions, branches,
                         branches.unresolved) {
  len <- cache$len
  depth <- cache$depth
  children <- cache$children
  order <- cache$order[-length(cache$order)]
  root <- cache$root

  unresolved <- cache$unresolved

  n <- length(len)
  lq <- rep(0, n)
  n.tip <- cache$n.tip

  y <- cache$y
  tips <- setdiff(seq_len(n.tip), unresolved$i)
  branch.init <- branch.base <- matrix(NA, n, ncol(y$y))

  if ( !is.null(unresolved) ) {
    i <- unresolved$i
    ans <- branches.unresolved(pars, len[i], unresolved)
    lq[i] <- ans[,1]
    branch.base[i,] <- ans[,-1]
  }

  if ( is.null(y$i) ) {
    branch.init[tips,] <- y$y
    for ( i in tips ) {
      ans <- branches(branch.init[i,], len[i], pars, 0)
      lq[i] <- ans[1]
      branch.base[i,] <- ans[-1]
    }
  } else {
    branch.init[tips,] <- y$y[y$i,]
    for ( i in y$types ) {
      idx <- tips[which(y$i == i)]
      t <- len[idx]
      ans <- branches(y$y[i,], sort(unique(t)), pars, 0)
      ans <- ans[tapply(t, t),,drop=FALSE]
      lq[idx] <- ans[,1]
      branch.base[idx,] <- ans[,-1]
    }
  }

  for ( i in order ) {
    y.in <- initial.conditions(branch.base[children[i,],], pars, depth[i])
    if ( !all(is.finite(y.in)) )
      stop("Bad initial conditions: calculation failure along branches?")
    branch.init[i,] <- y.in
    ans <- branches(y.in, len[i], pars, depth[i])
    lq[i] <- ans[1]
    branch.base[i,] <- ans[-1]
  }

  y.in <- initial.conditions(branch.base[children[root,],], pars,
                             TRUE, depth[root])
  branch.init[root,] <- y.in
  list(init=branch.init, base=branch.base, lq=lq)
}

## This is the minimal cache function, but not calculating initial
## conditions, which will take the element 'y'.
make.cache <- function(tree) {
  ## This works around some ape bugs with class inheritance.
  if (inherits(tree, "phylo"))
    class(tree) <- "phylo"
  edge <- tree$edge
  idx <- seq_len(max(edge))
  n.tip <- length(tree$tip.label)
  root <- n.tip + 1
  
  is.tip <- idx <= n.tip

  children <- lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
  if ( !all(sapply(children, length)==2) )
    stop("Multifircations/unbranched nodes in tree - best get rid of them")
  children <- rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), 2)))

  ## This is not carefully checked, but I think that it is right.
  parent <- rep(NA, length(idx))
  parent[-root] <- unlist(lapply(idx, function(x) edge[edge[,2]==x,1]))

  order <- get.ordering(children, is.tip, root)
  depth <- c(rep(0, n.tip), as.numeric(branching.times(tree)))

  anc <- ancestors(parent, order)
  
  ans <- list(len=tree$edge.len[match(idx, edge[,2])],
              children=children,
              parent=parent,
              order=order,
              root=root,
              n.tip=n.tip,
              depth=depth,
              ancestors=anc)
  ans
}

## Node ordering, as used by make.cache()
get.ordering <- function(children, is.tip, root) {
  todo <- list(root)
  i <- root
  repeat {
    kids <- children[i,]
    i <- kids[!is.tip[kids]]
    if ( length(i) > 0 )
      todo <- c(todo, list(i))
    else
      break
  }
  as.vector(unlist(rev(todo)))
}

ROOT.FLAT  <- 1
ROOT.EQUI  <- 2
ROOT.OBS   <- 3
ROOT.GIVEN <- 4
ROOT.BOTH  <- 5
ROOT.ALL   <- ROOT.BOTH
root.xxsse <- function(vars, pars, cache, condition.surv, root.mode,
                       root.p) {
  logcomp <- sum(vars$lq)
  vars <- vars$init[cache$root,]
  k <- length(vars) / 2
  e.root <- vars[seq_len(k)]
  d.root <- vars[(k+1):(2*k)]

  if ( root.mode == ROOT.FLAT )
    p <- 1/k
  else if ( root.mode == ROOT.OBS )
    p <- d.root / sum(d.root)
  else if ( root.mode == ROOT.EQUI )
    if ( k == 2 )
      p <- stationary.freq.bisse(pars)
    else
      stop("ROOT.EQUI only possible when k = 2")
  else if ( root.mode == ROOT.GIVEN ) {
    if ( length(root.p) != length(d.root) )
      stop("Invalid length for root.p")
    p <- root.p
  } else if ( root.mode != ROOT.ALL )
    stop("Invalid root mode")

  if ( condition.surv )
    d.root <- d.root / (1-e.root)^2
  if ( root.mode == ROOT.BOTH )
    loglik <- log(d.root)# + logcomp
  else
    loglik <- log(sum(p * d.root)) + logcomp
  loglik
}

cleanup <- function(loglik, pars, prior=NULL, intermediates=FALSE,
                    cache, vals) {
  if ( is.null(prior) )
    p <- loglik
  else if ( is.numeric(prior) )
    p <- loglik + prior.default(pars, prior)
  else if ( is.function(prior) )
    p <- loglik + prior(pars)
  else
    stop("Invalid 'prior' argument")
    
  if ( intermediates ) {
    attr(p, "cache") <- cache
    attr(p, "intermediates") <- vals
    attr(p, "vals") <- vals$init[cache$root,]
  }
  p
}

## Which leads to an all singing, all dancing function:
xxsse.ll <- function(pars, cache, initial.conditions,
                     branches, branches.unresolved, 
                     condition.surv, root.mode, root.p,
                     prior, intermediates) {
  ans <- all.branches(pars, cache, initial.conditions,
                      branches, branches.unresolved)
  loglik <- root.xxsse(ans, pars, cache, condition.surv,
                       root.mode, root.p)
  cleanup(loglik, pars, prior, intermediates, cache, ans)
}

make.prior.exponential <- function(r) {
  function(pars)
    -sum(pars * r)
}

prior.default <- function(pars, r) {
  .Deprecated("make.prior.exponential")
  - sum(pars * r)
}

## Convert a branches function into one that adds log-compensation.
## This is not compulsary to use, but should make life easier.
make.branches <- function(branches, idx, eps=0) {
  if ( length(idx) > 0 )
    function(y, len, pars, t0) {
      ret <- branches(y, len, pars, t0)
      if ( all(ret[,idx] >= eps) ) {
        q <- rowSums(ret[,idx,drop=FALSE])
        i <- q > 0
        ret[i,idx] <- ret[i,idx] / q[i]
        lq <- q
        lq[i] <- log(q[i])
        cbind(lq, ret, deparse.level=0)
      } else {
        ti <- len[length(len)]/2
        len1 <- c(len[len <= ti], ti)
        len2 <- len[len > ti] - ti
        n1 <- length(len1)
        ret1 <- Recall(y, len1, pars, t0)
        ret2 <- Recall(ret1[n1,-1], len2, pars, t0 + ti)
        ret2[,1] <- ret2[,1] + ret1[n1,1]
        rbind(ret1[-n1,], ret2)
      }
    }
  else
    function(y, len, pars, t0)
      cbind(0, branches(y, len, pars, t0), deparse.level=0)
}
