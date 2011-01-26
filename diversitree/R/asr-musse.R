## (1) Marginal ASR:  (untested)
asr.marginal.musse <- function(lik, pars, nodes=NULL, ...) {
  k <- attr(lik, "k") # I think, or possibly from cache?
  
  states.idx <- (k+1):(2*k)
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches
  root.p <- rep(1/k, k)  # Required
  condition.surv <- TRUE # Required (is this good or not?)

  res <- all.branches(pars, cache, initial.conditions.musse,
                      branches)

  root.f <- function(pars, vals, lq)
    root.xxsse(vals, pars, lq, condition.surv, root.p)
  
  do.asr.marginal(pars, cache, res, nodes, states.idx,
                  initial.conditions.musse,
                  branches, root.f)
}
