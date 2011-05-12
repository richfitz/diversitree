## (1) Marginal ASR:
asr.marginal.musse <- function(lik, pars, nodes=NULL,
                               condition.surv=TRUE,
                               root=ROOT.FLAT, root.p=NULL, ...) {
  k <- attr(lik, "k") # I think, or possibly from cache?
  states.idx <- (k+1):(2*k)
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches

  res <- all.branches(pars, cache, initial.conditions.musse,
                      branches)

  root.f <- function(pars, vals, lq)
    root.xxsse(vals, pars, lq, condition.surv,
               root.p.xxsse(vals, pars, root, root.p))
  
  do.asr.marginal(pars, cache, res, nodes, states.idx,
                  initial.conditions.musse,
                  branches, root.f)
}
