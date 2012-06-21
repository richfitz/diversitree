make.musse.t <- function(tree, states, k, functions, sampling.f=NULL,
                         strict=TRUE, control=list(),
                         spline.data=NULL) {
  cache <- make.cache.musse.t(tree, states, k, functions, sampling.f,
                              strict, spline.data)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.musse)
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
  class(ll) <- c("musse.t", "musse", "dtlik.t", "dtlik", "function")
  ll
}

make.cache.musse.t <- function(tree, states, k, functions,
                               sampling.f, strict, spline.data) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  update.cache.t(cache, functions, spline.data, with.q=TRUE)
}
