## (1) Marginal ASR:
asr.marginal.bisse <- function(lik, pars, nodes=NULL,
                               condition.surv=TRUE, root=ROOT.OBS,
                               root.p=NULL, ...) {
  states.idx <- 3:4
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches
  root.p <- rep(1/2, 2)

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
