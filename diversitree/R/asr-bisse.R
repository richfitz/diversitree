## (1) Marginal ASR:
asr.marginal.bisse <- function(lik, pars, nodes=NULL,
                               condition.surv=TRUE,
                               root=ROOT.FLAT, root.p=NULL, ...) {
  states.idx <- 3:4
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches

  if ( !is.null(cache$unresolved) )
    cache$preset <- 
      branches.unresolved.bisse(pars, cache$unresolved)

  res <- all.branches(pars, cache, initial.conditions.bisse,
                      branches)

  root.f <- function(pars, vals, lq)
    root.xxsse(vals, pars, lq, condition.surv,
               root.p.xxsse(vals, pars, root, root.p))
  
  do.asr.marginal(pars, cache, res, nodes, states.idx,
                  initial.conditions.bisse,
                  branches, root.f)
}
