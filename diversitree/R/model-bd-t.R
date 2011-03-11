## Time-varing birth-death models.  This is the simplest case, and
## will be good as a check, as I have the plain BD model implemented
## with the ODE system now.

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
make.bd.t <- function(tree, functions, sampling.f=NULL,
                      unresolved=NULL, control=list()) {
  control <- modifyList(list(safe=FALSE, tol=1e-8, eps=0), control)  
  cache <- make.cache.bd.ode(tree, sampling.f, unresolved)

  if ( is.null(names(functions)) && length(functions) == 2 )
    names(functions) <- argnames.bd(NULL)

  pars.t <- make.pars.t(functions)
  n.args <- attr(pars.t, "n.args")
  is.constant.arg <- attr(pars.t, "is.constant.arg")
  
  branches <- make.branches.bd.t(control$safe, control$tol,
                                 control$eps)
  initial.conditions <-
    make.initial.conditions.t(initial.conditions.bd.ode)
  const <- lfactorial(length(tree$tip.label) - 1)

  ll.bd.t <- function(pars, condition.surv=TRUE,
                      intermediates=FALSE) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid parameter length (expected %d)", n.args))
    pars.const <- pars[is.constant.arg]
    if ( any(pars.const < 0) || any(!is.finite(pars.const)) )
      return(-Inf)
    f.pars <- pars.t(pars)

    ll.xxsse.t(f.pars, cache, initial.conditions,
               branches, condition.surv, ROOT.FLAT, NULL,
               intermediates) + const
  }
  
  class(ll.bd.t) <- c("bd.t", "bd", "function")
  attr(ll.bd.t, "argnames") <- attr(pars.t, "argnames")
  ll.bd.t
}

`argnames<-.bd.t` <- function(x, value) {
  .NotYetImplemented()
}

make.branches.bd.t <- function(safe=FALSE, tol=1e-8, eps=0) {
  RTOL <- ATOL <- tol
  e <- new.env()

  bd.t <- make.ode("derivs_bd_t", "diversitree", "initmod_bd_t",
                       2, safe)
  branches <- function(y, len, pars, t0)
    t(bd.t(y, c(t0, t0+len), list(pars, e), rtol=RTOL,
           atol=ATOL)[-1,-1])

  make.branches(branches, 2, eps)
}

