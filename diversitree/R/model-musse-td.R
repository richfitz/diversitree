## 1: make
make.musse.td <- function(tree, states, k, n.epoch, sampling.f=NULL,
                          strict=TRUE, control=list()) {
  control <- check.control.ode(control)
  if ( control$backend == "CVODES" )
    stop("Cannot use CVODES backend with musse.td")

  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  cache$n.epoch <- n.epoch

  branches <- make.branches.td(make.branches.musse(cache, control))
  initial.conditions <-
    make.initial.conditions.td(initial.conditions.musse)

  npar <- (n.epoch - 1) + (k * (k + 1) * n.epoch)
  i.t <- seq_len(n.epoch - 1)
  i.p <- n.epoch:npar

  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.lm <- 1:(2*k)
  idx.q <- (2*k+1):(k*(1+k))

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    if ( length(pars) != npar )
      stop(sprintf("Invalid length parameters (expected %d)", npar))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    pars2 <- matrix(NA, n.epoch, k * (k + 2) + 1)
    pars2[,1] <- c(pars[i.t], Inf)
    tmp <- matrix(pars[i.p], n.epoch, k * (k + 1), TRUE)
    for ( i in seq_len(n.epoch) ) {
      qmat <- matrix(0, k, k)
      qmat[idx.qmat] <- tmp[i,idx.q]
      diag(qmat) <- -rowSums(qmat)
      pars2[i,-1] <- c(tmp[i,idx.lm], qmat)
    }

    ll.xxsse.td(pars2, cache, initial.conditions, branches,
                condition.surv, root, root.p, intermediates)
  }
  
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
