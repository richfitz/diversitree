make.bd.t.old <- function(tree, functions, sampling.f=NULL,  
                          unresolved=NULL, control=list()) {
  if (interactive()) warning("make.bd.t.old is temporary only")
  cache <- make.cache.bd.t.old(tree, functions, sampling.f, unresolved)
  all.branches <-
    make.all.branches.t.old.dtlik(cache, control,
                                  initial.conditions.bd.ode)
  rootfunc <- make.rootfunc.t.old(cache, rootfunc.bd.ode)
  pars.t <- make.pars.t.old(cache$functions, cache)
  const <- cache$const

  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    f.pars <- pars.t(pars)
    ans <- all.branches(f.pars, intermediates)
    rootfunc(ans, f.pars, condition.surv, intermediates, const)
  }
  class(ll) <- c("bd.t", "bd", "dtlik", "function")
  ll
}

make.cache.bd.t.old <- function(tree, functions, unresolved, sampling.f) { 
  cache <- make.cache.bd.ode(tree, unresolved, sampling.f)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  cache$info$ml.default <- "subplex"
  update.cache.t.old(cache, functions)
  
}

## Not used directly.
make.branches.bd.t.old <- function(cache, control)
  make.branches.dtlik.t.old(cache$info, control)
