## (1) Marginal ASR:
make.asr.marginal.musse.split <- function(lik, ...) {
  e <- environment(lik)
  k <- e$k
  cache <- e$cache
  branches <- e$branches
  initial.conditions <- e$initial.conditions
  f.pars <- e$f.pars
  n.part <- cache$n.part
  unresolved <- cache$unresolved

  use.CVODES <- is.null(branches)

  states.idx <- (k+1):(2*k)

  if ( use.CVODES )
    stop("asr.marginal.musse.split with CVODES not yet implemented")

  g.root <- cache$group.nodes[cache$root] # always 1 but being safe.

  function(pars, nodes=NULL, condition.surv=TRUE,
           root=ROOT.FLAT, root.p=NULL, ...) {
    root.f <- function(pars, vals, lq)
      root.xxsse(vals, pars[[g.root]], lq, condition.surv,
                 root.p.xxsse(vals, pars[[g.root]], root, root.p))

    pars <- check.par.musse.split(pars, n.part, k)
    for ( i in seq_len(n.part) )
      pars[[i]] <- f.pars(pars[[i]])

    if ( use.CVODES ) {
      stop("Unpossible!")
    } else {
      res <- all.branches.matrix(pars, cache, initial.conditions, branches)
      do.asr.marginal(pars, cache, res, nodes, states.idx,
                      initial.conditions,
                      branches, root.f)
    }
  }
}
