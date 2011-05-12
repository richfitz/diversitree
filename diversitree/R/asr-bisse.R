## (1) Marginal ASR:
asr.marginal.bisse <- function(lik, pars, nodes=NULL,
                               condition.surv=TRUE,
                               root=ROOT.FLAT, root.p=NULL, ...) {
  e <- environment(lik)
  states.idx <- 3:4
  cache <- e$cache
  branches <- e$branches

  if ( is.null(branches) ) {
    control <- e$control
    if ( control$backend == "CVODES" ) {
      control$backend <- "cvodes"
      branches <- make.branches.bisse(cache, control)
    } else {
      stop("'branches' missing from likelihood function")
    }
  }

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
