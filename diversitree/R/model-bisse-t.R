
make.bisse.t <- function(tree, states, functions, sampling.f=NULL,
                         unresolved=NULL, nt.extra=10, strict=TRUE,
                         safe=FALSE) {
  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra,
                            strict=strict)

  if ( is.null(names(functions)) && length(functions) == 6 )
    names(functions) <- argnames.bisse(NULL)
  pars.t <- make.pars.t(functions)
  n.args <- attr(pars.t, "n.args")
  is.constant.arg <- attr(pars.t, "is.constant.arg")

  branches <- make.branches.bisse.t(safe)
  initial.conditions <-
    make.initial.conditions.t(initial.conditions.bisse)

  ll.bisse.t <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                         root.p=NULL, intermediates=FALSE) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid length parameters (expected %d)", n.args))
    pars.const <- pars[is.constant.arg]
    if ( any(pars.const < 0) || any(!is.finite(pars.const)) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    f.pars <- pars.t(pars)

    ll.xxsse.t(f.pars, cache, initial.conditions, branches,
               condition.surv, root, root.p, intermediates)
  }

  class(ll.bisse.t) <- c("bisse.t", "bisse", "function")
  attr(ll.bisse.t, "argnames") <- attr(pars.t, "argnames")
  ll.bisse.t
}

`argnames<-.bisse.t` <- function(x, value) {
  .NotYetImplemented()
}

## 8: branches
make.branches.bisse.t <- function(safe=FALSE) {
  RTOL <- ATOL <- 1e-8
  e <- new.env()
  
  bisse.t.ode <- make.ode("derivs_bisse_t", "diversitree",
                          "initmod_bisse_t", 4, safe)
  branches <- function(y, len, pars, t0)
    t(bisse.t.ode(y, c(t0, t0+len), list(pars, e),
                  rtol=RTOL, atol=ATOL)[-1,-1])
  make.branches(branches, 3:4)
}
