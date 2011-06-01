make.asr.marginal.mkn <- function(lik, ...) {
  k <- as.integer(attr(lik, "k"))
  states.idx <- seq_len(k)
  cache <- environment(lik)$cache
  len <- cache$len

  cache.C <- list(parent=toC.int(cache$parent),
                  children=toC.int(t(cache$children)),
                  root=toC.int(cache$root))
  env <- new.env()
  nodes.all.C <- toC.int(cache$root:max(cache$order))

  function(pars, nodes=NULL, root=ROOT.FLAT, root.p=NULL, ...) {
    root.f <- function(pars, vals, lq)
      root.mkn(vals, lq,
               root.p.mkn(vals, pars, root, root.p))

    qmat <- mkn.Q(pars)
    res <- all.branches.mkn(qmat, cache)

    ## TODO: ugly!  This is basically duplicated in do.asr.marginal()
    ## and do.asr.marginal.C(); there should be nicer way of doing
    ## this.
    if ( is.null(nodes) )
      nodes.C <- nodes.all.C
    else
      nodes.C <- toC.int(nodes + cache$n.tip)

    .Call("r_asr_marginal_mkn", k, pars, nodes.C, cache.C, res,
          root.f, env, PACKAGE="diversitree")
  }
}

make.asr.joint.mkn <- function(lik, ...) {
  k <- as.integer(attr(lik, "k"))
  cache <- environment(lik)$cache
  is.mk2 <- inherits(lik, "mk2")

  order.C <- toC.int(rev(cache$order))
  parent.C <- toC.int(cache$parent)

  function(pars, n=1, ...) {
    obj <- attr(lik(pars, intermediates=TRUE, ...), "intermediates")
    ##.Call("r_do_asr_joint", k, order.C, parent.C, obj$init, obj$pij,
    ##      obj$root.p, as.01, PACKAGE="diversitree")
    do.asr.joint(n, k, order.C, parent.C, obj$init, obj$pij,
                 obj$root.p, is.mk2)
  }
}

make.asr.stoch.mkn <- function(lik, slim=FALSE, ...) {
  is.mk2 <- inherits(lik, "mk2")
  cache <- environment(lik)$cache
  k <- as.integer(attr(lik, "k"))

  if ( is.mk2 ) {
    tip.state <- as.integer(cache$tip.state - 1L)
    pos.states <- c(0L, 1L)
  } else {
    tip.state <- as.integer(cache$tip.state)
    pos.states <- seq_len(k)
  }

  edge.1 <- cache$edge[,1]
  edge.2 <- cache$edge[,2]
  edge.length <- cache$edge.length # == cache$len[edge.2]

  ## Joint distribution
  joint <- make.asr.joint(lik)

  ## Single branch simulator
  ptr <- .Call("r_smkn_alloc", k, 100L, PACKAGE="diversitree")
  sim.set.pars <- function(pars)
    .Call("r_smkn_set_pars", ptr, pars, PACKAGE="diversitree")
  sim1 <- function(len, state.beg, state.end)
    .Call("r_smkn_scm_run", ptr, len, state.beg, state.end, is.mk2,
          PACKAGE="diversitree")

  ## TODO: This is fragile, and should be stored in the cache anyway.
  n.node <- nrow(cache$edge)/2

  function(pars, n=1, ...) {
    if ( n > 1 )
      stop("Not yet implemented (n>1)")

    ## 1: Draw a random sample from the nodes' joint distribution:
    node.state <- joint(pars, n, ...)

    ## If we are using mkn, we need to deflate the node states onto
    ## base-0 indices for the simulation code.
    if ( is.mk2 )
      anc.state <- as.integer(c(tip.state, node.state))
    else
      anc.state <- as.integer(c(tip.state - 1L, node.state - 1L))

    ## 2: Simulate branches.
    ## a: set the parameters; the same values will be used for all
    ## branches.
    sim.set.pars(pars)

    ## b: Sort the state beginning, end and extinction into the order
    ## that we will use.
    state.beg <- anc.state[edge.1]
    state.end <- anc.state[edge.2]

    ## c: Actually simulate the histories
    f <- function(i)
      sim1(edge.length[i], state.beg[i], state.end[i])
    history <- lapply(seq_along(state.beg), f)

    if ( slim ) {
      tmp <- lapply(history, function(x) x[-1,,drop=FALSE])
      keep <- which(sapply(tmp, length) > 0)
      names(node.state) <- NULL
      ret <- list(node.state=node.state,
                  history=list(idx=keep, tmp[keep]))
    } else {
      if ( !is.null(cache$node.label) )
        names(node.state) <- cache$node.label
      make.history(NULL, tip.state, node.state, history, TRUE,
                   pos.states, check=FALSE)
    }
  }
}

summarise.histories.mk2 <- function(x, phy) {
  summarise.branch <- function(x) {
    n <- length(x)
    j <- seq_len(n)
    y <- cbind(do.call(rbind, x),
               rep(j, sapply(x, nrow)), deparse.level=0)
    y <- y[order(y[,1]),]
    cbind(c(0, y[-j,1]),
          (c(0, cumsum(y[-j,2] * 2 - 1)) + sum(y[j,2])) / n,
          deparse.level=0)
  }

  nn <- length(x[[1]]$history)
  nx <- length(x)
  tmp <- matrix(unlist(lapply(x, "[[", "history"), FALSE), nn, nx)
  h <- apply(tmp, 1, summarise.branch)

  node.state <- rowMeans(matrix(unlist(lapply(x, "[[", "node.state")),
                                nn/2, nx))
  names(node.state) <- names(x[[1]]$node.state)
  tip.state <- x[[1]]$tip.state
  make.history(phy, tip.state, node.state, h, FALSE, 0:1)
}
