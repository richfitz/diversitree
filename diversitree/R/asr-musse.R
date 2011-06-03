## (1) Marginal ASR:
make.asr.marginal.musse <- function(lik, ...) {
  e <- environment(lik)  
  k <- attr(lik, "k")
  states.idx <- (k+1):(2*k)
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches
  use.CVODES <- is.null(branches)

  if ( use.CVODES ) {
    control <- e$control
    if ( control$backend != "CVODES" )
      stop("'branches' missing from likelihood function")      
    all.branches.C <- e$all.branches
    ptr <- environment(all.branches.C)$ptr
    env <- new.env()
    states.idx.C <- toC.int(states.idx)
    parent.C <- toC.int(cache$parent)
  }

  f.pars <- make.musse.pars(k)

  function(pars, nodes=NULL, condition.surv=TRUE,
           root=ROOT.FLAT, root.p=NULL, ...) {
    root.f <- function(pars, vals, lq)
      root.xxsse(vals, pars, lq, condition.surv,
                 root.p.xxsse(vals, pars, root, root.p))
    check.pars.musse(pars, k)
    pars2 <- f.pars(pars)

    if ( use.CVODES ) {
      do.asr.marginal.C(pars2, cache, ptr, nodes, states.idx.C,
                        parent.C, all.branches.C, root.f, env)
    } else {
      res <- all.branches.matrix(pars2, cache,
                                 initial.conditions.musse,
                                 branches)
      do.asr.marginal(pars2, cache, res, nodes, states.idx,
                      initial.conditions.musse,
                      branches, root.f)
    }
  }
}
