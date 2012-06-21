make.bd.t <- function(tree, functions, sampling.f=NULL,  
                      unresolved=NULL, control=list(),
                      spline.data=NULL) {
  cache <- make.cache.bd.t(tree, functions, sampling.f, unresolved,
                           spline.data)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.bd.ode)
  rootfunc <- make.rootfunc.t(cache, rootfunc.bd.ode)
  const <- cache$const

  ## TODO: This is currently uglier than necessary.  
  prep <-
    make.prep.all.branches.t(cache,
                             environment(all.branches)$control$backend)

  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    ans <- all.branches(prep(pars), intermediates)
    rootfunc(ans, pars, condition.surv, intermediates, const)
  }
  class(ll) <- c("bd.t", "bd", "dtlik.t", "dtlik", "function")
  ll
}

make.cache.bd.t <- function(tree, functions, unresolved, sampling.f,
                            spline.data) {
  cache <- make.cache.bd.ode(tree, unresolved, sampling.f)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  cache$info$ml.default <- "subplex"
  update.cache.t(cache, functions, spline.data)
}
