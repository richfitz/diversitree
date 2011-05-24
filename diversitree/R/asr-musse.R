## (1) Marginal ASR:
asr.marginal.musse <- function(lik, ...) {
  make.asr.marginal(lik)(...)
}

make.asr.marginal.musse <- function(lik, ...) {
  e <- environment(lik)  
  k <- attr(lik, "k") # I think, or possibly from cache?
  states.idx <- (k+1):(2*k)
  cache <- environment(lik)$cache
  branches <- environment(lik)$branches

  if ( is.null(branches) ) {
    control <- e$control
    if ( control$backend != "CVODES" ) {
      stop("'branches' missing from likelihood function")      
    } else if ( use.CVODES ) {
      all.branches.C <- e$all.branches
      ptr <- environment(all.branches.C)$ptr
      env <- new.env()
      states.idx.C <- toC.int(states.idx)
      parent.C <- toC.int(cache$parent)
    } else {
      control$backend <- "cvodes"
      branches <- make.branches.musse(cache, control)
    }
  } else {
    use.CVODES <- FALSE
  }

  function(pars, nodes=NULL,
           condition.surv=TRUE,
           root=ROOT.FLAT, root.p=NULL, ...) {
    root.f <- function(pars, vals, lq)
      root.xxsse(vals, pars, lq, condition.surv,
                 root.p.xxsse(vals, pars, root, root.p))

    if ( use.CVODES ) {
      do.asr.marginal.C(pars, cache, ptr, nodes, states.idx.C,
                        parent.C, all.branches.C, root.f, env)
    } else {
      res <- all.branches(pars, cache, initial.conditions.musse,
                          branches)
      do.asr.marginal(pars, cache, res, nodes, states.idx,
                      initial.conditions.musse,
                      branches, root.f)
    }
  }
}
