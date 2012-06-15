make.bd.t2 <- function(tree, functions, sampling.f=NULL,  
                       unresolved=NULL, control=list(),
                       spline.data=NULL) {
  cache <- make.cache.bd.t2(tree, functions, sampling.f, unresolved,
                            spline.data)
  all.branches <- make.all.branches.t2.dtlik(cache, control,
                                            initial.conditions.bd.ode)
  rootfunc <- make.rootfunc.t2(cache, rootfunc.bd.ode)

  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, intermediates)
  }
  class(ll) <- c("bd.t2", "bd", "dtlik", "function")
  ll
}

make.cache.bd.t2 <- function(tree, functions, unresolved, sampling.f,
                             spline.data) {
  cache <- make.cache.bd.ode(tree, unresolved, sampling.f)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  cache$info$ml.default <- "subplex"
  update.cache.t2(cache, functions, spline.data)
}
