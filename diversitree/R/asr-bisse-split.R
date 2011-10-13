## (1) Marginal ASR:
make.asr.marginal.bisse.split <- function(lik, ...) {
  e <- environment(lik)
  states.idx <- 3:4
  cache <- e$cache
  branches <- e$branches
  initial.conditions <- e$initial.conditions
  n.part <- cache$n.part
  unresolved <- cache$unresolved

  use.CVODES <- is.null(branches)

  if ( use.CVODES )
    stop("asr.marginal.bisse.split with CVODES not yet implemented")

  g.root <- cache$group.nodes[cache$root] # always 1 but being safe.

  function(pars, nodes=NULL, condition.surv=TRUE,
           root=ROOT.FLAT, root.p=NULL, ...) {
    root.f <- function(pars, vals, lq)
      root.xxsse(vals, pars[[g.root]], lq, condition.surv,
                 root.p.xxsse(vals, pars[[g.root]], root, root.p))

    pars <- check.par.bisse.split(pars, n.part)

    if ( use.CVODES ) {
      stop("Impossible!")
    } else {
      if ( !is.null(unresolved) )
        cache$preset <- 
          branches.unresolved.bisse.split(pars, unresolved)

      res <- all.branches.matrix(pars, cache, initial.conditions, branches)
      do.asr.marginal(pars, cache, res, nodes, states.idx,
                      initial.conditions,
                      branches, root.f)
    }
  }
}
