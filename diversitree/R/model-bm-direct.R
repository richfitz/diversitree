## Simple-minded BM calculator, direct from the solution to the BM
## PDF.
## 1: make
make.bm.direct <- function(tree, states, states.sd=0) {
  cache <- make.cache.bm.direct(tree, states, states.sd)

  ll.bm.direct <- function(pars, root=ROOT.MAX, root.x=NA,
                     intermediates=FALSE) {
    check.pars.bm(pars)
    ans <- all.branches.matrix(pars, cache, initial.conditions.bm.direct,
                               branches.bm.direct)
    vals <- ans$init[,cache$root]
    loglik <- root.bm.direct(vals, ans$lq, root, root.x)
    if ( intermediates ) {
      attr(loglik, "intermediates") <- intermediates
      attr(loglik, "vals") <- vals
    }
    loglik
  }
  
  class(ll.bm.direct) <- c("bm.direct", "bm", "function")
  ll.bm.direct
}

## 2: print (inherited from bm)
## 3: argnames / argnames<- (inherited from bm)
argnames.bm.direct <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    "s2"
  else
    ret
}

## 4: find.mle (inherited from bm)

## 5: make.cache
make.cache.bm.direct <- function(tree, states, states.sd) {
  tree <- check.tree(tree, ultrametric=FALSE)

  if ( is.null(states.sd) )
    states.sd <- 0
  tmp <- check.states.quasse(tree, states, states.sd)
  states <- tmp$states
  states.sd <- tmp$states.sd

  cache <- make.cache(tree)
  cache$ny <- 3
  y <- mapply(function(mean, sd) c(mean, sd*sd, 0),
              states, states.sd, SIMPLIFY=TRUE)
  cache$y.tmp <- y
  cache$y <- dt.tips.ordered(y, cache$tips, cache$len[cache$tips])

  cache
}

## 6: ll
root.bm.direct <- function(vars, log.comp, root, root.x) {
  if ( root == ROOT.MAX ) {
    ## This treats a prior on the root as a delta function centred at
    ## the ML root state.
    ## The first term can be more intuitively written as:
    ##   dnorm(vars[1], vars[1], sqrt(vars[2]), TRUE)
    ##   dnorm(0, 0, sqrt(vars[2]), TRUE)
    - log(2 * pi * vars[[2]]) / 2 + vars[[3]] + sum(log.comp)
  } else if ( root == ROOT.FLAT ) {
    ## Flat prior (by this point, function integrates to vars[[3]])
    vars[[3]] + sum(log.comp)
  } else if ( root == ROOT.OBS ) {
    ## Observed weighting (integrate norm norm wrt x from -inf to inf
    ## gives 1 / (2 sqrt(pi s2))).
    -log(2 * sqrt(pi * vars[[2]])) + vars[[3]] + sum(log.comp)
  } else if ( root == ROOT.GIVEN ) {
    dnorm(root.x, vars[1], sqrt(vars[2]), TRUE) + vars[[3]] +
      sum(log.comp)
  } else {
    stop("Invalid root mode")
  }
}

## 7: initial.conditions
initial.conditions.bm.direct <- function(init, pars, t, idx) {
  m1 <- init[1,1]
  m2 <- init[1,2]
  v1 <- init[2,1]
  v2 <- init[2,2]
  vv <- v1 + v2

  c((m1 * v2 + m2 * v1) / vv,
    (v1 * v2) / vv,
    -(m1 - m2)^2 / (2 * vv) - log(2 * pi * vv) / 2)
}

## 8: branches
branches.bm.direct <- function(y, len, pars, t0, idx) {
  m <- y[1]
  v <- y[2]
  z <- y[3]
  list(z, c(m, v + (pars[1] * len), 0))
}
