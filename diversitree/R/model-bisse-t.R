
make.bisse.t <- function(tree, states, functions, sampling.f=NULL,
                         unresolved=NULL, nt.extra=10, strict=TRUE,
                         control=list()) {
  control <- check.control.ode(control)
  if ( control$backend == "CVODES" )
    stop("Cannot use CVODES backend with bisse.t")

  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra,
                            strict=strict)

  if ( is.null(names(functions)) && length(functions) == 6 )
    names(functions) <- argnames.bisse(NULL)
  pars.t <- make.pars.t(functions)
  n.args <- attr(pars.t, "n.args")
  is.constant.arg <- attr(pars.t, "is.constant.arg")

  branches <- make.branches.bisse.t(cache, control)
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
make.branches.bisse.t <- function(cache, control) {
  neq <- 4L
  np <- 6L # not used
  comp.idx <- as.integer(3:4)
  make.ode.branches.t("bisse_t", "diversitree", neq, np, comp.idx,
                      control)
}
