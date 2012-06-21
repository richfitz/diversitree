make.bisse.t <- function(tree, states, functions, unresolved=NULL,
                         sampling.f=NULL, strict=TRUE, control=list(),
                         spline.data=NULL) {
  cache <- make.cache.bisse.t(tree, states, functions, unresolved,
                              sampling.f, strict, spline.data)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.bisse)
  rootfunc <- make.rootfunc.t(cache, rootfunc.musse)

  ## TODO: This is currently uglier than necessary.  
  prep <-
    make.prep.all.branches.t(cache,
                             environment(all.branches)$control$backend)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    ans <- all.branches(prep(pars), intermediates)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.t", "bisse", "dtlik.t", "dtlik", "function")
  ll
}

make.cache.bisse.t <- function(tree, states, functions,
                               unresolved, sampling.f,
                               strict, spline.data) {
  nt.extra <- 10 # fixed, as we're going to ignore unresolved.
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  update.cache.t(cache, functions, spline.data)
}
