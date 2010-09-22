## 1: make
make.bisse.td <- function(tree, states, n.epoch, unresolved=NULL,
                          sampling.f=NULL, nt.extra=10, strict=TRUE,
                          safe=FALSE) {
  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra,
                            strict=strict)
  cache$n.epoch <- n.epoch
  if ( !is.null(cache$unresolved) )
    stop("Cannot (yet) use unresolved clades with time-dependent BiSSE")

  branches <- make.branches.td(make.branches.bisse(safe))
  initial.conditions <-
    make.initial.conditions.td(initial.conditions.bisse)

  npar <- (n.epoch - 1) + (6 * n.epoch)
  i.t <- seq_len(n.epoch - 1)
  i.p <- n.epoch:npar

  ll.bisse.td <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {

    if ( length(pars) != npar )
      stop(sprintf("Invalid length parameters (expected %d)", npar))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    pars <- cbind(c(pars[i.t], Inf),
                  matrix(pars[i.p], n.epoch, 6, TRUE))

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
  c(sprintf("t.%d", seq_len(n.epoch-1)),
    argnames.twopart(x, argnames.bisse(NULL), n.epoch))
}
`argnames<-.bisse.td` <- function(x, value) {
  n.epoch <- attr(x, "n.epoch")
  argnames.twopart.set(x, value, 6, n.epoch)
}

## 4: find.mle: from bisse

## Make requires the usual functions:
## 5: make.cache (in model-bisse)

## 6: ll

## 7: initial.conditions

## 8: branches
