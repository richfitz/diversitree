## This is the general calculation code for diversitree.  Many models
## require that calculations occur down a branch, and the results of
## these calculations are used as initial conditions for the next
## branch.  There are several different implemented models using this
## approach:
##   - BiSSE (model-bisse.R)
##   - MuSSE (multi-state BiSSE - model-musse.R)
##   - QuaSSE
##   - Time-dependent BiSSE, MuSSE
##   - BD2 (HIV model)

## I have abstracted most of the calculations in a (hopefully) general
## way.  The function 'all.branches' does most of the calculations;
## this looks after all the book-keeping as calculations proceed down
## a tree.  Models will generally use this as their main function, but
## with additional root state calculations after this calculation has
## finished (see model-bisse.R for the canonical example of this).

## all.branches takes arguments
##   pars, cache, initial.conditions, branches
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
all.branches <- function(pars, cache, initial.conditions, branches) {
  ## Inside all.branches:
  len <- cache$len
  depth <- cache$depth
  children <- cache$children
  order <- cache$order[-length(cache$order)]
  root <- cache$root

  n <- length(len)
  lq <- rep(0, n)
  n.tip <- cache$n.tip

  y <- cache$y
  branch.init <- branch.base <- vector("list", n)

  ## TODO: This should also move in the tip conditions perhaps?
  if ( !is.null(cache$preset) ) {
    lq[cache$preset$target] <- cache$preset$lq
    branch.base[cache$preset$target] <- cache$preset$base
  }

  ## TODO: better way of sorting between these.
  ## how about a tip.method=TIP.SORTED or TIP.GROUPED?
  if ( is.null(names(y)) ) {
    for ( x in y ) {
      idx <- x$target
      branch.init[idx] <- list(x$y)
      ## TODO: the 'zero' here assumes that the tree is ultrametric.
      ## I can pass in
      ##   cache$depth[idx]
      ## but I need to be careful with this, as some depths will be
      ## 1e-15 and things like that.
      ans <- branches(x$y, x$t.uniq, pars, 0)

      if ( is.matrix(ans) ) {
        ans <- ans[x$unpack,,drop=FALSE]
        lq[idx] <- ans[,1]
        branch.base[idx] <- matrix.to.list(ans[,-1,drop=FALSE])
      } else {
        ans <- ans[x$unpack]
        lq[idx] <- unlist(lapply(ans, "[[", 1))
        branch.base[idx] <- lapply(ans, "[", -1)
      }
    }
  } else {
    tip.y <- branch.init[cache$tips] <- y$y
    tip.t <- y$t
    tip.target <- y$target
    for ( i in seq_along(tip.y) ) {
      ## TODO: This (the zero) assumes that all branches terminate at
      ## the present (not true for Mk2 style models, or extinct
      ## species)
      j <- tip.target[i]
      ans <- branches(tip.y[[i]], tip.t[i], pars, 0)
      lq[j] <- ans[1]
      branch.base[[j]] <- ans[-1]
    }    
  }

  for ( i in order ) {
    y.in <- initial.conditions(branch.base[children[i,]], pars, depth[i])
    if ( !all(is.finite(y.in)) )
      stop("Bad initial conditions: calculation failure along branches?")
    branch.init[[i]] <- y.in
    ans <- branches(y.in, len[i], pars, depth[i])
    lq[i] <- ans[1]
    branch.base[[i]] <- ans[-1]
  }

  ## This also changes to reflect the change in argument order.
  y.in <- initial.conditions(branch.base[children[root,]], pars,
                             depth[root], TRUE)
  branch.init[[root]] <- y.in
  list(init=branch.init, base=branch.base, lq=lq)
}

## This is the minimal cache function, but not calculating initial
## conditions, which will take the element 'y'.
make.cache <- function(tree) {
  ## This works around some ape bugs with class inheritance.
  if (inherits(tree, "phylo"))
    class(tree) <- "phylo"
  edge <- tree$edge
  edge.length <- tree$edge.length
  idx <- seq_len(max(edge))
  n.tip <- length(tree$tip.label)
  tips <- seq_len(n.tip)
  root <- n.tip + 1
  
  is.tip <- idx <= n.tip

  children <- lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
  ## Technically, this is already checked by check.tree, but I'm happy
  ## leaving it in.
  if ( !all(sapply(children, length)==2) )
    stop("Multifircations/unbranched nodes in tree - best get rid of them")
  children <- rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), 2)))

  parent <- edge[match(idx, edge[,2]),1]

  order <- get.ordering(children, is.tip, root)
  len <- edge.length[match(idx, edge[,2])]

  ## This is a bit of a hack, but this is to ensure that we can really
  ## compute the depths accurately - this is a problem when there
  ## joins (under split models) that occur right around nodes.
  height <- branching.heights(tree)
  depth <- max(height) - height
  depth2 <- branching.depth(len, children, order, tips)
  i <- abs(depth - depth2) < 1e-8
  depth[i] <- depth2[i]

  ## TODO: I don't need this ancestor thing for much - drop it here
  ## and move it to the asr code that actually uses it (this takes a
  ## lot of time, and is only used by the ASR code).
  anc <- vector("list", max(order))
  for ( i in c(rev(order[-length(order)]), tips) )
    anc[[i]] <- c(parent[i], anc[[parent[i]]])
  
  ans <- list(tip.label=tree$tip.label,
              len=len,
              children=children,
              parent=parent,
              order=order,
              root=root,
              n.tip=n.tip,
              tips=tips,
              height=height,
              depth=depth,
              ancestors=anc,
              edge=edge,
              edge.length=edge.length)
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

root.p.xxsse <- function(vals, pars, root, root.p=NULL) {
  k <- length(vals) / 2
  d.root <- vals[(k+1):(2*k)]

  if ( root == ROOT.FLAT )
    p <- 1/k
  else if ( root == ROOT.EQUI )
    if ( k == 2 )
      p <- stationary.freq.bisse(pars)
    else
      stop("ROOT.EQUI only possible when k = 2")
  else if ( root == ROOT.OBS )
    p <- d.root / sum(d.root)
  else if ( root == ROOT.GIVEN ) {
    if ( length(root.p) != length(d.root) )
      stop("Invalid length for root.p")
    p <- root.p
  } else if ( root == ROOT.ALL )
    p <- NULL
  else
    stop("Invalid root mode")
  p
}

root.xxsse <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)

  k <- length(vals) / 2
  i <- seq_len(k)
  lambda <- pars[i]
  e.root <- vals[i]
  d.root <- vals[-i]
  
  if ( condition.surv )
    d.root <- d.root / (lambda * (1-e.root)^2)

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

cleanup <- function(loglik, pars, intermediates=FALSE, cache, vals) {
  if ( intermediates ) {
    attr(loglik, "cache") <- cache
    attr(loglik, "intermediates") <- vals
    attr(loglik, "vals") <- vals$init[[cache$root]]
  }
  loglik
}

## Which leads to an all singing, all dancing function:
ll.xxsse <- function(pars, cache, initial.conditions,
                     branches, condition.surv, root, root.p,
                     intermediates) {
  ans <- all.branches(pars, cache, initial.conditions, branches)
  vals <- ans$init[[cache$root]]
  root.p <- root.p.xxsse(vals, pars, root, root.p)
  loglik <- root.xxsse(vals, pars, ans$lq, condition.surv, root.p)
  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
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

## Utility functions for organising initial conditions.
## TODO: Document.
dt.tips.grouped <- function(y, y.i, tips, t) {
  if ( !is.list(y) )
    stop("'y' must be a list of initial conditions")
  if ( max(y.i) > length(y) || min(y.i) < 1 )
    stop("'y.i' must be integers on 1..", length(y))
  if ( length(y.i) != length(tips) )
    stop("y must be same length as tips")
  if ( length(y.i) != length(t) )
    stop("y must be the same length as t")
  
  types <- sort(unique(y.i))
  res <- vector("list", length(types))

  for ( i in seq_along(types) ) {
    j <- which(y.i == types[i])
    target <- tips[j]
    t.i <- t[j]
    t.uniq <- sort(unique(t.i))
    unpack <- match(t.i, t.uniq)
    res[[i]] <- list(y=y[[types[i]]], y.i=types[i], target=target,
                     t.uniq=t.uniq, unpack=unpack)
  }
  res
}

dt.tips.ordered <- function(y, tips, t) {
  if ( !is.list(y) )
    stop("'y' must be a list of initial conditions")
  if ( length(y) != length(tips) )
    stop("y must be same length as tips")
  if ( length(y) != length(t) )
    stop("y must be the same length as t")

  i <- order(t)
  list(target=tips[i],
       t=t[i],
       y=y[i])
}
