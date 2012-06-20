make.bisse.t.old <- function(tree, states, functions, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10, strict=TRUE,
                             control=list()) {
  if (interactive()) warning("make.bisse.t.old is temporary only")  
  cache <- make.cache.bisse.t.old(tree, states, functions, unresolved,
                                  sampling.f, nt.extra, strict)
  all.branches <- make.all.branches.t.old.dtlik(cache, control,
                                                initial.conditions.bisse)
  rootfunc <- make.rootfunc.t.old(cache, rootfunc.musse)
  pars.t <- make.pars.t.old(cache$functions, cache)
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    f.pars <- pars.t(pars)
    ans <- all.branches(f.pars, intermediates)
    rootfunc(ans, f.pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.t", "bisse", "dtlik", "function")
  ll
}

make.cache.bisse.t.old <- function(tree, states, functions,
                                   unresolved, sampling.f,
                                   nt.extra, strict) {
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying BiSSE with unresolved clades")
  update.cache.t.old(cache, functions)
}

## Not used directly.
make.branches.bisse.t.old <- function(cache, control)
  make.branches.dtlik.t.old(cache$info, control)
