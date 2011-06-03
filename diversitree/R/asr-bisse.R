## (1) Marginal ASR:
make.asr.marginal.bisse <- function(lik, ...) {
  e <- environment(lik)
  states.idx <- 3:4
  cache <- e$cache
  branches <- e$branches

  use.CVODES <- is.null(branches)

  if ( is.null(branches) ) {
    control <- e$control
    if ( control$backend != "CVODES" )
      stop("'branches' missing from likelihood function")      

    all.branches.C <- e$all.branches
    ptr <- environment(all.branches.C)$ptr
    env <- new.env()
    states.idx.C <- toC.int(states.idx)
    parent.C <- toC.int(cache$parent)
  }

  function(pars, nodes=NULL, condition.surv=TRUE,
           root=ROOT.FLAT, root.p=NULL, ...) {
    root.f <- function(pars, vals, lq)
      root.xxsse(vals, pars, lq, condition.surv,
                 root.p.xxsse(vals, pars, root, root.p))
    if ( !check.pars.bisse(pars) )
      stop("Invalid parameters")

    if ( use.CVODES ) {
      do.asr.marginal.C(pars, cache, ptr, nodes, states.idx.C,
                        parent.C, all.branches.C, root.f, env)
    } else {
      if ( !is.null(cache$unresolved) )
        cache$preset <- 
          branches.unresolved.bisse(pars, cache$unresolved)

      res <- all.branches.matrix(pars, cache,
                                 initial.conditions.bisse,
                                 branches)
      do.asr.marginal(pars, cache, res, nodes, states.idx,
                      initial.conditions.bisse,
                      branches, root.f)
    }
  }
}
