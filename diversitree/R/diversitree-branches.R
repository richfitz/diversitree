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

## The 'as.list' argument is to handle cases, such as QuaSSE, where
## the output from branches() is of variable length and the result
## should be stored in a list, rather than a matrix.  It is *not* to
## switch between a choice of what is returned.

## TODO: the 'zero' in the tip branches() calls assume that the tree
## is ultrametric.  I can pass in
##   cache$depth[idx]
## but I need to be careful with this, as some depths will be
## 1e-15 and things like that.  This may not be a problem in reality.
all.branches.matrix <- function(pars, cache, initial.conditions,
                                branches) {
  len <- cache$len
  depth <- cache$depth
  children <- cache$children
  order <- cache$order[-length(cache$order)]
  root <- cache$root

  n <- length(len)
  lq <- rep(0, n)
  n.tip <- cache$n.tip

  y <- cache$y
  branch.init <- branch.base <- matrix(NA, cache$ny, n)

  if ( !is.null(cache$preset) ) {
    lq[cache$preset$target] <- cache$preset$lq
    branch.base[,cache$preset$target] <- cache$preset$base
  }

  if ( is.null(names(y)) ) { # dt.tips.grouped
    for ( x in y ) {
      if ( !is.null(x) ) {
        ## The above is because sometimes 'x' will be NULL.  This
        ## happens when there is no polymorphism across tips, usually
        ## when most of the tips are unresolved clades.
        idx <- x$target
        branch.init[,idx] <- x$y
        ans <- branches(x$y, x$t, pars, 0, idx)
        lq[idx] <- ans[[1]]
        branch.base[,idx] <- ans[[2]]
      }
    }
  } else { # y$type == "ORDERED"
    tip.t <- y$t
    tip.target <- y$target
    tip.y <- branch.init[,tip.target] <- y$y
    for ( i in seq_along(tip.t) ) {
      idx <- tip.target[i]
      ans <- branches(tip.y[,i], tip.t[i], pars, 0, idx)
      lq[idx] <- ans[[1]]
      branch.base[,idx] <- ans[[2]]
    }
  }

  for ( i in order ) {
    y.in <- initial.conditions(branch.base[,children[i,]], pars,
                               depth[i], i)
    if ( !all(is.finite(y.in)) )
      stop("Bad initial conditions: calculation failure along branches?")
    branch.init[,i] <- y.in
    ans <- branches(y.in, len[i], pars, depth[i], i)
    lq[i] <- ans[[1]]
    branch.base[,i] <- ans[[2]]
  }

  y.in <- initial.conditions(branch.base[,children[root,]], pars,
                             depth[root], root)
  branch.init[,root] <- y.in
  list(init=branch.init, base=branch.base, lq=lq)
}

all.branches.list <- function(pars, cache, initial.conditions,
                              branches) {
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

  if ( !is.null(cache$preset) ) {
    lq[cache$preset$target] <- cache$preset$lq
    branch.base[cache$preset$target] <- cache$preset$base
  }

  if ( is.null(names(y)) ) { # dt.tips.grouped
    for ( x in y ) {
      if ( !is.null(x) ) {
        idx <- x$target
        branch.init[idx] <- list(x$y)
        ans <- branches(x$y, x$t, pars, 0, idx)
        lq[idx] <- unlist(lapply(ans, "[[", 1))
        branch.base[idx] <- lapply(ans, "[", -1)
      }
    }
  } else { # dt.tips.ordered
    tip.t <- y$t
    tip.target <- y$target
    tip.y <- branch.init[tip.target] <- y$y
    for ( i in seq_along(tip.t) ) {
      idx <- tip.target[i]
      ans <- branches(tip.y[[i]], tip.t[i], pars, 0, idx)
      lq[idx] <- ans[1]
      branch.base[[idx]] <- ans[-1]
    }
  }

  for ( i in order ) {
    y.in <- initial.conditions(branch.base[children[i,]], pars,
                               depth[i], i)
    if ( !all(is.finite(y.in)) )
      stop("Bad initial conditions: calculation failure along branches?")
    branch.init[[i]] <- y.in
    ans <- branches(y.in, len[i], pars, depth[i], i)
    lq[i] <- ans[1]
    branch.base[[i]] <- ans[-1]
  }

  y.in <- initial.conditions(branch.base[children[root,]], pars,
                             depth[root], root)
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
  ## The only place that this is used at all is do.asr.marginal(); it
  ## would be possible to make this as needed when making an
  ## asr.marginal() function.
  anc <- vector("list", max(order))
  for ( i in c(rev(order[-length(order)]), tips) )
    anc[[i]] <- c(parent[i], anc[[parent[i]]])

  ans <- list(tip.label=tree$tip.label,
              node.label=tree$node.label,
              len=len,
              children=children,
              parent=parent,
              order=order,
              root=root,
              n.tip=n.tip,
              n.node=tree$Nnode,
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
ROOT.MAX   <- 6
ROOT.ALL   <- ROOT.BOTH

root.p.xxsse <- function(vals, pars, root, root.p=NULL) {
  k <- length(vals) / 2
  d.root <- vals[(k+1):(2*k)]

  if ( root == ROOT.FLAT ) {
    p <- 1/k
  } else if ( root == ROOT.EQUI ) {
    if ( k == 2 ) {
      ## TODO: This is a bit of an ugliness now that other models have
      ## stationary frequencies (bisseness, geosse, classe).
      p <- stationary.freq.bisse(pars)
      p <- c(p, 1-p)
    } else {
      stop("ROOT.EQUI only implemented for BiSSE")
    }
  } else if ( root == ROOT.OBS ) {
    p <- d.root / sum(d.root)
  } else if ( root == ROOT.GIVEN ) {
    if ( length(root.p) != length(d.root) )
      stop("Invalid length for root.p")
    p <- root.p
  } else if ( root == ROOT.ALL ) {
    p <- NULL
  } else {
    stop("Invalid root mode")
  }
  p
}

root.xxsse <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)
  is.root.both <- is.null(root.p)

  k <- length(vals) / 2
  i <- seq_len(k)
  lambda <- pars[i]
  e.root <- vals[i]
  d.root <- vals[-i]

  if ( condition.surv ) {
    if ( is.root.both )
      d.root <- d.root / (lambda * (1-e.root)^2)
    else
      d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2)
  }

  if ( is.root.both ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

## Which leads to an all singing, all dancing function:
ll.xxsse <- function(pars, cache, initial.conditions,
                     branches, condition.surv, root, root.p,
                     intermediates) {
  ans <- all.branches.matrix(pars, cache, initial.conditions, branches)
  vals <- ans$init[,cache$root]
  root.p <- root.p.xxsse(vals, pars, root, root.p)
  loglik <- root.xxsse(vals, pars, ans$lq, condition.surv, root.p)

  if ( intermediates ) {
    ans$root.p <- root.p
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- vals
  }

  loglik
}

## Convert a branches function into one that adds log-compensation.
## This is not compulsary to use, but should make life easier.
make.branches <- function(branches, comp.idx, eps=0) {
  if ( length(comp.idx) > 0 )
    function(y, len, pars, t0, idx) {
      ret <- branches(y, len, pars, t0, idx)
      q <- colSums(ret[comp.idx,,drop=FALSE])
      if ( all(q >= eps) ) {
        i <- q > 0
        ret[comp.idx,i] <- ret[comp.idx,i] /
          rep(q[i], each=length(comp.idx))
        lq <- q
        lq[i] <- log(q[i])
        list(lq, ret)
      } else {
        ti <- len[length(len)]/2
        len1 <- c(len[len <= ti], ti)
        len2 <- len[len > ti] - ti
        n1 <- length(len1)

        ret1 <- Recall(y, len1, pars, t0)
        ret2 <- Recall(ret1[[2]][,n1], len2, pars, t0 + ti)
        ret2[[1]] <- ret2[[1]] + ret1[[1]][n1]

        list(c(ret1[[1]][-n1], ret2[[1]]),
             cbind(ret1[[2]][,-n1], ret2[[2]]))
      }
    }
  else
    function(y, len, pars, t0, idx)
      list(rep.int(0, length(len)),
           branches(y, len, pars, t0, idx))
}

make.ode.branches <- function(model, dll, neq, np, comp.idx, control) {
  backend <- control$backend
  safe <- control$safe
  tol <- control$tol
  eps <- control$eps

  if ( backend == "deSolve" ) {
    initfunc <- sprintf("initmod_%s", model)
    derivs <- sprintf("derivs_%s", model)

    RTOL <- ATOL <- tol
    ode <- make.ode(derivs, dll, initfunc, neq, safe)
    branches <- function(y, len, pars, t0, idx)
      ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1,drop=FALSE]
  } else if ( backend == "cvodes" ) {
    derivs <- sprintf("derivs_%s_cvode", model)

    RTOL <- ATOL <- tol
    ode <- cvodes(neq, np, derivs, RTOL, ATOL, dll)
    branches <- function(y, len, pars, t0, idx)
      ode(pars, y, c(t0, t0+len))[,-1,drop=FALSE]
  } else {
    stop("Invalid backend", backend)
  }

  make.branches(branches, comp.idx, eps)
}

## Utility functions for organising initial conditions.
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

  for ( type in types ) {
    j <- which(y.i == type)
    i <- order(t[j])
    res[[type]] <- list(y=y[[type]], y.i=type,
                        target=tips[j][i], t=t[j][i])
  }
  res
}

dt.tips.ordered <- function(y, tips, t) {
  i <- order(t)

  if ( is.list(y) ) {
    if ( length(y) != length(tips) )
      stop("y must be same length as tips")
    if ( length(y) != length(t) )
      stop("y must be the same length as t")
    list(target=tips[i],
         t=t[i],
         y=y[i])
  } else if ( is.matrix(y) ) {
    if ( ncol(y) != length(tips) )
      stop("y must be same length as tips")
    if ( ncol(y) != length(t) )
      stop("y must be the same length as t")
    list(target=tips[i],
         t=t[i],
         y=y[,i], type="ORDERED")
  } else {
    stop("y must be a list or matrix")
  }
}

