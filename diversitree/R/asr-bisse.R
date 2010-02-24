## (1) Marginal ASR:
asr.marginal.bisse <- function(lik, pars, nodes=NULL, ...) {
  states.idx <- 3:4
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches
  root.p <- rep(1/2, 2)
  condition.surv <- TRUE # TODO: fix pending...

  res <- all.branches(pars, cache, initial.conditions.bisse,
                      branches, branches.unresolved.bisse)

  root.f <- function(pars, vals, lq)
    root.xxsse(vals, pars, lq, condition.surv, root.p)
  
  do.asr.marginal(pars, cache, res, nodes, states.idx,
                  initial.conditions.bisse,
                  branches,
                  branches.unresolved.bisse,
                  root.f)
}
