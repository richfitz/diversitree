## (1) Marginal ASR:
asr.marginal.bisse <- function(lik, pars, nodes=NULL, ...) {
  states.idx <- 3:4
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches

  res <- all.branches(pars, cache, initial.conditions.bisse,
                      branches, branches.unresolved.bisse)
  
  do.asr.marginal(pars, cache, res, nodes, states.idx,
                  initial.conditions.bisse,
                  branches,
                  branches.unresolved.bisse)
}
