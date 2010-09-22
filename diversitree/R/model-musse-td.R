## 1: make
make.musse.td <- function(tree, states, k, n.epoch, sampling.f=NULL,
                          strict=TRUE, safe=FALSE) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  cache$n.epoch <- n.epoch

  branches <- make.branches.td(make.branches.musse(k, safe))

  initial.conditions <-
    make.initial.conditions.td(initial.conditions.musse)

  npar <- (n.epoch - 1) + (k * (k + 1) * n.epoch)
  i.t <- seq_len(n.epoch - 1)
  i.p <- n.epoch:npar

  ll.musse.td <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                          root.p=NULL, intermediates=FALSE) {
    if ( length(pars) != npar )
      stop(sprintf("Invalid length parameters (expected %d)", npar))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    pars <- cbind(c(pars[i.t], Inf),
                  matrix(pars[i.p], n.epoch, k * (k + 1), TRUE))

    ll.xxsse.td(pars, cache, initial.conditions, branches,
                condition.surv, root, root.p, intermediates)
  }
  
  ll <- function(pars, ...) ll.musse.td(pars, ...)

  class(ll) <- c("musse.td", "musse", "function")
  attr(ll, "n.epoch") <- n.epoch
  attr(ll, "k") <- k
  ll
}

## 2: print
print.musse.td <- function(x, ...) {
  cat("MuSSE/td likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.musse.td <- function(x, k=attr(x, "k"),
                                 n.epoch=attr(x, "n.epoch"), ...) {
  c(sprintf("t.%d", seq_len(n.epoch-1)),
    argnames.twopart(x, argnames.musse(NULL, k), n.epoch))
}
`argnames<-.musse.td` <- function(x, value) {
  n.epoch <- attr(x, "n.epoch")
  k <- attr(x, "k")
  argnames.twopart.set(x, value, k * (k + 1), n.epoch)
}

## 4: find.mle
find.mle.musse.td <- function(func, x.init, method, fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.musse.td")
}

## Make requires the usual functions:
## 5: make.cache (in model-musse)

## 6: ll

## 7: initial.conditions

## 8: branches
