make.bd.t <- function(tree, functions, sampling.f=NULL,  
                      unresolved=NULL, control=list()) {
  cache <- make.cache.bd.t(tree, functions, sampling.f, unresolved)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.bd.ode)
  rootfunc <- make.rootfunc.t(cache, rootfunc.bd.ode)
  pars.t <- make.pars.t(cache$functions, cache)

  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    f.pars <- pars.t(pars)
    ans <- all.branches(f.pars, intermediates)
    rootfunc(ans, f.pars, condition.surv, intermediates)
  }
  class(ll) <- c("bd.t", "bd", "dtlik", "function")
  ll
}

make.cache.bd.t <- function(tree, functions, unresolved, sampling.f) { 
  cache <- make.cache.bd.ode(tree, unresolved, sampling.f)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  cache$info$ml.default <- "subplex"
  update.cache.t(cache, functions)
  
}

## Not used directly.
make.branches.bd.t <- function(cache, control)
  make.branches.dtlik.t(cache$info, control)
