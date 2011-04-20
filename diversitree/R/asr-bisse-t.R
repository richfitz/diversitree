asr.marginal.bisse.t <- function(lik, pars, nodes=NULL, ...) {
  states.idx <- 3:4
  e <- environment(lik)
  
  cache <- e$cache
  branches <- e$branches
  branches <- e$branches
  initial.conditions <- e$initial.conditions
  
  root.p <- rep(1/2, 2)  # TODO: Thinking required
  condition.surv <- TRUE # TODO: fix pending...

  if ( !is.null(cache$unresolved) )
    stop("Unresolved/time varying bisse/asr marginal not yet done")

  f.pars <- e$pars.t(pars)

  pars.root <- f.pars(cache$depth[cache$root])
  root.f <- function(pars, vals, lq)
    root.xxsse(vals, pars.root, lq, condition.surv, root.p)

  res <- all.branches(f.pars, cache, initial.conditions,
                      branches)
  
  do.asr.marginal(f.pars, cache, res, nodes, states.idx,
                  initial.conditions,
                  branches, root.f)
}

asr.marginal.bisse.td <- function(lik, pars, nodes=NULL, ...) {
  stop("Not yet implemented.")
}
