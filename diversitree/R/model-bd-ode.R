## ODE interface to BD models.

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
make.bd.ode <- function(tree, sampling.f=NULL, unresolved=NULL,
                        control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend

  cache <- make.cache.bd.ode(tree, sampling.f, unresolved)

  if ( backend == "CVODES" )
    all.branches <- make.all.branches.C.bd(cache, control)
  else
    branches <- make.branches.bd(cache, control)
    
  const <- lfactorial(length(tree$tip.label) - 1)

  ll.bd.ode <- function(pars, condition.surv=TRUE,
                        intermediates=FALSE) {
    if ( length(pars) != 2 )
      stop("Invalid parameter length (expected 2)")
    if ( any(pars < 0) || any(!is.finite(pars)) )
      return(-Inf)

    if ( backend == "CVODES" )
      ll.xxsse.C(pars, all.branches,
                 condition.surv, ROOT.FLAT, NULL,
                 intermediates) + const
    else
      ll.xxsse(pars, cache, initial.conditions.bd.ode,
               branches, condition.surv, ROOT.FLAT, NULL,
               intermediates) + const
  }
  
  class(ll.bd.ode) <- c("bd.ode", "bd", "function")
  ll.bd.ode
}

## 2: print: from bd
## 3: argnames / argnames<-: from bd
## 4: find.mle: from bd

## 5: make.cache (initial.tip, root)
make.cache.bd.ode <- function(tree, unresolved, sampling.f) {
  tree <- check.tree(tree)

  if ( inherits(tree, "clade.tree") || !is.null(unresolved) )
    stop("Cannot deal with unresolved clades yet with this method.")

  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, 1)

  cache <- make.cache(tree)
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f
  cache$y <- initial.tip.bd.ode(cache)
  cache
}

initial.tip.bd.ode <- function(cache) {
  f <- cache$sampling.f
  tips <- cache$tips
  
  y <- list(c(1-f, f))
  y.i <- rep(1, length(tips))
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}

## 6: ll: internal

## 7: initial.conditions:
initial.conditions.bd.ode <- function(init, pars, t, is.root=FALSE)
  c(init[[1]][1], init[[1]][2] * init[[2]][2] * pars[1])

## 8: branches
make.branches.bd <- function(cache, control) {
  neq <- 2L
  np <- 2L
  comp.idx <- 2L
  make.ode.branches("bd", "diversitree", neq, np, comp.idx, control)
}

make.all.branches.C.bd <- function(cache, control) {
  neq <- 2L
  np <- 2L
  comp.idx <- 2L
  make.all.branches.C(cache, "bd", "diversitree",
                      neq, np, comp.idx, control)
}
