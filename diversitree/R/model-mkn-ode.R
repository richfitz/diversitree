## ODE interface to MK models.

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
make.mkn.ode <- function(tree, states, k, strict=TRUE, control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend

  cache <- make.cache.mkn.ode(tree, states, k, strict)

  if ( backend == "CVODES" )
    all.branches <- make.all.branches.C.mkn.ode(cache, control)
  else
    branches <- make.branches.mkn.ode(cache, control)

  f.pars <- make.mkn.pars(k)  
    
  ll <- function(pars, root=ROOT.OBS, root.p=NULL,
                 intermediates=FALSE) {
    if (!is.null(root.p) && root != ROOT.GIVEN) 
      warning("Ignoring specified root state")
    check.pars.mkn(pars, k)

    qmat <- f.pars(pars)
    if (backend == "CVODES") {
      ans <- all.branches(qmat, intermediates)
      vals <- ans[[2]]
      root.p <- root.p.mkn(vals, pars, root, root.p)
      loglik <- root.mkn(vals, ans[[1]], root.p)
      if ( intermediates ) {
        ans$intermediates$root.p <- root.p
        attr(loglik, "intermediates") <- ans$intermediates
        attr(loglik, "vals") <- vals
      }
    } else {
      ans <- all.branches.matrix(qmat, cache,
                                 initial.conditions.mkn.ode, branches)
      vals <- ans$init[,cache$root]
      root.p <- root.p.mkn(vals, pars, root, root.p)
      loglik <- root.mkn(vals, ans$lq, root.p)
      if ( intermediates ) {
        ans$root.p <- root.p
        attr(loglik, "intermediates") <- ans
        attr(loglik, "vals") <- vals
      }
    }
    loglik
  }
  
  class(ll) <- c("mkn.ode", "mkn", "function")
  attr(ll, "k") <- k    
  ll
}

## 2: print: from mkn
## 3: argnames / argnames<-: from mkn
## 4: find.mle: from mkn

## 5: make.cache (initial.tip, root)
make.cache.mkn.ode <- function(tree, states, k, strict=TRUE) {
  tree <- check.tree(tree)

  states <- check.states(tree, states,
                         strict=strict, strict.vals=1:k)
  states <- check.integer(states)

  cache <- make.cache(tree)
  cache$ny <- k
  cache$k <- k
  cache$tip.state <- states
  cache$y <- initial.tip.mkn.ode(cache)
  cache
}

initial.tip.mkn.ode <- function(cache) {
  k <- cache$k

  y <- matrix(rep(c(rep(0, k)), k + 1), k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,]) <- 1
  y <- matrix.to.list(y)

  y.i <- cache$tip.state
  y.i[is.na(y.i)] <- k + 1

  tips <- cache$tips
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}

## 6: ll: internal

## 7: initial.conditions:
initial.conditions.mkn.ode <- function(init, pars, t, idx) {
  ret <- init[,1] * init[,2]
  if ( !any(ret > 0) )
    stop("Incompatible initial conditions at tip ", idx)
  ret
}

## 8: branches
make.branches.mkn.ode <- function(cache, control) {
  k <- as.integer(cache$k)
  neq <- k
  np <- k*k
  comp.idx <- seq_len(k)
  make.ode.branches("mknode", "diversitree", neq, np, comp.idx, control)
}

make.all.branches.C.mkn.ode <- function(cache, control) {
  k <- as.integer(cache$k)
  neq <- k
  np <- k*k
  comp.idx <- seq_len(k)
  make.all.branches.C(cache, "mknode", "diversitree",
                      neq, np, comp.idx, control)
}
