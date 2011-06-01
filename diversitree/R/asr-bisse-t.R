make.asr.marginal.bisse.t <- function(lik, ...) {
  e <- environment(lik)
  states.idx <- 3:4
  cache <- e$cache
  branches <- e$branches
  initial.conditions <- e$initial.conditions

  if ( is.null(branches) )
    ## Should not even get here, but it future proofs us.
    stop("Not possible with CVODES backend")

  function(pars, nodes=NULL, condition.surv=TRUE,
           root=ROOT.FLAT, root.p=NULL, ...) {
    f.pars <- e$pars.t(pars)
    pars.root <- f.pars(cache$depth[cache$root])
    
    root.f <- function(pars, vals, lq)
      root.xxsse(vals, pars, lq, condition.surv,
                 root.p.xxsse(vals, pars, root, root.p))

    res <- all.branches.matrix(f.pars, cache, initial.conditions,
                               branches)
    do.asr.marginal(f.pars, cache, res, nodes, states.idx,
                    initial.conditions,
                    branches, root.f)
  }
}

asr.marginal.bisse.td <- function(lik, pars, nodes=NULL, ...) {
  stop("Not yet implemented.")
}
