make.bisse.t2 <- function(tree, states, functions, unresolved=NULL,
                          sampling.f=NULL, strict=TRUE, control=list(),
                          spline.data=NULL) {
  cache <- make.cache.bisse.t2(tree, states, functions, unresolved,
                               sampling.f, strict, spline.data)
  all.branches <- make.all.branches.t2.dtlik(cache, control,
                                             initial.conditions.bisse)
  rootfunc <- make.rootfunc.t2(cache, rootfunc.musse)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.t2", "bisse", "dtlik.t", "dtlik", "function")
  ll
}

make.cache.bisse.t2 <- function(tree, states, functions,
                                unresolved, sampling.f,
                                strict, spline.data) {
  nt.extra <- 10 # fixed, as we're going to ignore unresolved.
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  update.cache.t2(cache, functions, spline.data)
}
