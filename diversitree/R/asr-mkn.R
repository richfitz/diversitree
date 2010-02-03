asr.marginal.mkn <- function(lik, pars, nodes=NULL, ...) {
  k <- attr(lik, "k")
  states.idx <- seq_len(k)
  cache <- environment(lik)$cache
  len <- cache$len

  qmat <- mkn.Q(pars)
  res <- all.branches.mkn(qmat, cache)
  pij <- res$pij

  branches <- function(y, len.i, pars, t0) {
    if ( length(len.i) != 1 )
      stop("Should not happen.")
    res <- t(matrix(pij[,match(len.i, len)], k, k) %*% y)
    q <- rowSums(res)
    cbind(log(q), res/q, deparse.level=0)
  }
  
  do.asr.marginal(pars, cache, res, nodes, states.idx,
                  initial.conditions.mkn,
                  branches,
                  branches.unresolved.mkn)
}

asr.joint.mkn <- function(lik, pars, n=1, root.state=NA,
                          simplify=TRUE, intermediates=FALSE,
                          ...) {
  ## Current:
  k <- attr(lik, "k")
  cache <- environment(lik)$cache
  node.labels <- environment(lik)$tree$node.label
  obj <- attr(lik(pars, intermediates=TRUE), "intermediates")
  li <- obj$init
  pij <- t(obj$pij)

  ## Newer/better?
  ##   k <- attr(lik, "k")
  ##   cache <- environment(lik)$cache
  ##   qmat <- mkn.Q(pars)
  ##   obj <- all.branches.mkn(qmat, cache)
  ##   li <- obj$init  
  ##   pij <- t(obj$pij)

  is.mk2 <- inherits(lik, "mk2")
  if ( is.mk2 && !is.na(root.state) )
    root.state <- root.state - 1

  x <- do.asr.joint(pars, n, root.state, cache, li, pij,
                    node.labels, simplify)
  
  if ( is.mk2 ) {
    if ( !simplify )
      x[] <- lapply(x, function(el) el - 1)
    else
      x <- x - 1
  }

  if ( intermediates )
    attr(x, "intermediates") <- obj
  x
}

asr.stoch.mkn <- function(lik, pars, n=1, root.state=NA, ...) {
  is.mk2 <- inherits(lik, "mk2")
  k <- attr(lik, "k")
  cache <- environment(lik)$cache
  edge <- cache$edge
  edge.length <- cache$edge.length

  node.state <- asr.joint(lik, pars, n, root.state,
                          simplify=TRUE, intermediates=TRUE)

  if ( n == 1 )
    do.asr.stoch.mkn.one(pars, cache$tip.state, node.state,
                         edge, edge.length, k, is.mk2)
  else
    replicate(n,
              do.asr.stoch.mkn.one(pars, cache$tip.state, node.state,
                                   edge, edge.length, k, is.mk2),
              simplify=FALSE)
}

do.asr.stoch.mkn.one <- function(pars, tip.state, node.state,
                                 edge, edge.length, k, as.01) {
  if ( as.01 ) {
    anc.state <- c(tip.state, node.state + 1)
    tip.state <- tip.state - 1
    states <- c(0, 1)
  } else {
    anc.state <- c(tip.state, node.state)
    states <- 1:k
  }
  
  state.beg <- anc.state[edge[,1]]
  state.end <- anc.state[edge[,2]]

  f <- function(i)
    stoch.branch.mkn(pars, edge.length[i],
                     state.beg[i], state.end[i], k, as.01)

  history <- lapply(seq_along(state.beg), f)
  make.history(NULL, tip.state, node.state, history, TRUE, states,
               FALSE)
}

stoch.branch.mkn <- function(pars, len, state.beg, state.end, k,
                             as.01) {
  pars <- matrix(pars, k, k-1, TRUE)
  q.diag <- rowSums(pars) 
  idx <- seq_len(k)
  i <- lapply(idx, function(i) idx[-i])

  ## TODO: There is also a special case for three states, but I bet that
  ## this can be done faster in C for any number of states.  Probably
  ## also changing the first time to sample from a conditional
  ## exponential when state.beg != state.end?
  f2 <- function(state) {
    t <- 0
    changes <- list(c(t, state))
    while ( (t <- t + rexp(1, pars[state])) < len ) {
      state <- i[[state]]
      changes[[length(changes)+1]] <- c(t, state)
    }
    list(state, changes)
  }
  
  f3 <- function(state) {
    t <- 0
    changes <- list(c(t, state))
    while ( (t <- t + rexp(1, pars[state])) < len ) {
      state <- i[[state]][sample(k-1, 1, FALSE, pars[1,])]
      changes[[length(changes)+1]] <- c(t, state)      
    }
    list(state, changes)    
  }

  f <- if ( k == 2 ) f2 else f3

  repeat {
    tmp <- f(state.beg)
    if ( state.end == tmp[[1]] )
      break
  }

  ans <- do.call(rbind, tmp[[2]])
  if ( as.01 )
    ans[,2] <- ans[,2] - 1
  ans
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
