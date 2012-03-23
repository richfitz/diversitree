make.bisse.t <- function(tree, states, functions, unresolved=NULL,
                         sampling.f=NULL, nt.extra=10, strict=TRUE,
                         control=list()) {
  cache <- make.cache.bisse.t(tree, states, functions, unresolved,
                              sampling.f, nt.extra, strict)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.bisse)
  rootfunc <- make.rootfunc.t(cache, rootfunc.musse)
  pars.t <- make.pars.t(cache$functions, cache)
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    f.pars <- pars.t(pars)
    ans <- all.branches(f.pars, intermediates)
    rootfunc(ans, f.pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.t", "bisse", "dtlik", "function")
  ll
}

make.cache.bisse.t <- function(tree, states, functions,
                               unresolved, sampling.f,
                               nt.extra, strict) {
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying BiSSE with unresolved clades")
  update.cache.t(cache, functions)
}

## Not used directly.
make.branches.bisse.t <- function(cache, control)
  make.branches.dtlik.t(cache$info, control)
