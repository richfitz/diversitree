## 1: make
make.bisse.td <- function(tree, states, n.epoch, unresolved=NULL,
                          sampling.f=NULL, nt.extra=10, safe=FALSE) {
  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra)
  cache$n.epoch <- n.epoch
  if ( !is.null(cache$unresolved) )
    stop("Cannot (yet) use unresolved clades with time-dependent BiSSE")

  branches <- make.branches.td(make.branches.bisse(safe))
  initial.conditions <-
    make.initial.conditions.td(initial.conditions.bisse)

  ll.bisse.td <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    n.epoch <- cache$n.epoch
    if ( length(pars) != n.epoch * 7 )
      stop("Expected pars of length ", 7 * n.epoch)
    pars <- matrix(pars, n.epoch, 7, TRUE)
    if ( any(pars < 0) || any(!is.finite(pars[,-1])) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    ll.xxsse.td(pars, cache, initial.conditions, branches,
                condition.surv, root, root.p, intermediates)
  }
  
  ll <- function(pars, ...) ll.bisse.td(pars, ...)
  class(ll) <- c("bisse.td", "bisse", "function")
  attr(ll, "n.epoch") <- n.epoch
  ll
}

## 2: print
print.bisse.td <- function(x, ...) {
  cat("BiSSE/td likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.bisse.td <- function(x, n.epoch=attr(x, "n.epoch"), ...) {
  obj <- attr(x, "argnames")
  if ( is.null(obj) ) {
    obj <- list(par=c("t", "lambda0", "lambda1", "mu0", "mu1", "q01", "q10"),
                t=sprintf("t%d", 1:n.epoch))
  }
  paste(obj$par, rep(obj$t, each=7), sep=".")
}
`argnames<-.bisse.td` <- function(x, value) {
  .NotYetImplemented()
}

## 4: find.mle
find.mle.bisse.td <- function(func, x.init, method, fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.bisse.td")
}

## Make requires the usual functions:
## 5: make.cache (in model-bisse)

## 6: ll

## 7: initial.conditions

## 8: branches
