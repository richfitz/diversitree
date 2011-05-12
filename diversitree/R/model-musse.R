## Models should provide:
##   1. make
##   2. print
##   3. argnames / argnames<-
##   4. find.mle
## Generally, make will require:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches

## 1: make
make.musse <- function(tree, states, k, sampling.f=NULL, strict=TRUE,
                       control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend

  cache <- make.cache.musse(tree, states, k, sampling.f, strict)

  if ( backend == "CVODES" )
    all.branches <- make.all.branches.C.musse(cache, control)
  else
    branches <- make.branches.musse(cache, control)

  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.lm <- 1:(2*k)
  idx.q <- (2*k+1):(k*(1+k))
  
  ll.musse <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                       root.p=NULL, intermediates=FALSE) {
    if ( length(pars) != k*(k+1) )
      stop(sprintf("Invalid length parameters (expected %d)",
                   k*(k+1)))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    pars2 <- c(pars[idx.lm], qmat)

    if ( backend == "CVODES" )
      ll.xxsse.C(pars2, all.branches,
                 condition.surv, root, root.p, intermediates)
    else
      ll.xxsse(pars2, cache, initial.conditions.musse, branches,
               condition.surv, root, root.p, intermediates)
  }

  ll <- function(pars, ...) ll.musse(pars, ...)
  class(ll) <- c("musse", "function")
  attr(ll, "k") <- k
  ll
}

## 2: print
print.musse <- function(x, ...) {
  cat("MuSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.musse <- function(x, k=attr(x, "k"), ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    fmt <- sprintf("%%0%dd", ceiling(log10(k + .5)))
    str <- sprintf(fmt, 1:k)
    c(sprintf("lambda%s", str),
      sprintf("mu%s", str),
      sprintf("q%s%s", rep(str, each=k-1),
              unlist(lapply(1:k, function(i) str[-i]))))
  } else {
    ret
  }
}
`argnames<-.musse` <- function(x, value) {
  k <- attr(x, "k")
  ## k <- environment(x)$cache$k  
  if ( length(value) != k*(k+1) )
    stop("Invalid names length")
  if ( any(duplicated(value)) )
    stop("Duplicate argument names")
  attr(x, "argnames") <- value
  x
}

## 4: find.mle
find.mle.musse <- function(func, x.init, method, fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.musse")
}

mcmc.musse <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache (initial.tip)
make.cache.musse <- function(tree, states, k, sampling.f=NULL,
                             strict=TRUE) {
  tree <- check.tree(tree)
  states <- check.states(tree, states,
                         strict=strict, strict.vals=1:k)
  states <- check.integer(states)

  sampling.f <- check.sampling.f(sampling.f, k)

  cache <- make.cache(tree)
  cache$tip.state <- states
  cache$k <- k
  cache$sampling.f <- sampling.f
  cache$y <- initial.tip.musse(cache)
  cache
}
initial.tip.musse <- function(cache) {
  k <- cache$k
  f <- cache$sampling.f

  y <- matrix(rep(c(1-f, rep(0, k)), k + 1), k+1, 2*k, TRUE)
  y[k+1,(k+1):(2*k)] <- diag(y[1:k,(k+1):(2*k)]) <- f
  y <- matrix.to.list(y)

  y.i <- cache$tip.state
  y.i[is.na(y.i)] <- k + 1

  if ( !is.null(multistate <- attr(cache$tip.state, "multistate")) ) {
    y.multi <- unique(multistate$states)
    y.i.multi <- match(multistate$states, y.multi)

    y <- c(y, lapply(y.multi, function(x) c(1-f, x)))
    y.i[multistate$i] <- y.i.multi + k + 1
  }

  tips <- cache$tips
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}

## 6: ll.musse is done within make.musse

## 7: initial.conditions:
initial.conditions.musse <- function(init, pars, t, is.root=FALSE) {
  k <- length(init[[1]])/2
  i <- seq_len(k)
  j <- i + k

  c(init[[1]][i],
    init[[1]][j] * init[[2]][j] * pars[i])
}

## 8: branches
make.branches.musse <- function(cache, control) {
  k <- cache$k
  if ( k == 2 )
    warning("Two states is faster with BiSSE")
  np <- as.integer(k * (k + 2))
  neq <- as.integer(2 * k)
  comp.idx <- as.integer((k+1):(2*k))
  make.ode.branches("musse", "diversitree", neq, np, comp.idx,
                    control)
}

make.all.branches.C.musse <- function(cache, control) {
  k <- cache$k
  np <- as.integer(k * (k + 2))
  neq <- as.integer(2 * k)
  comp.idx <- as.integer((k+1):(2*k))
 
  make.all.branches.C(cache, "musse", "diversitree",
                      neq, np, comp.idx, control)
}

## Historical interest: This function creates a function for computing
## derivatives under MuSSE.  However, it uses the wrong argument setup
## right now, so I've commented it out.  It probably should not have
## been used.
## make.musse.eqs.R <- function(k) {
##   qmat <- matrix(0, k, k)
##   idx.qmat <- cbind(rep(1:k, each=k-1),
##                unlist(lapply(1:k, function(i) (1:k)[-i])))
##   idx.e <- 1:k
##   idx.d <- (k+1):(2*k)
##   idx.l <- 1:k
##   idx.m <- (k+1):(2*k)
##   idx.q <- (2*k+1):(k*(1+k))
##   function(t, y, parms, ...) {
##     e <- y[idx.e]
##     d <- y[idx.d]
##     lambda <- parms[idx.l]
##     mu     <- parms[idx.m]
##     qmat[idx.qmat] <- parms[idx.q]  
##     diag(qmat) <- -rowSums(qmat)
##     list(c(mu - (lambda + mu) * e + lambda * e * e + qmat %*% e,
##            -(lambda + mu) * d + 2 * lambda * d * e + qmat %*% d))
##   }
## }

## This makes the Q matrix from a set of parameters.
musse.Q <- function(pars, k) {
  if ( missing(k) )
    k <- (sqrt(1+4*length(pars))-1)/2
  if ( abs(k - round(k)) > 1e-6 || length(pars) != k*(1+k) )
    stop("Invalid parameter length")
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.q <- (2*k+1):(k*(1+k))
  qmat[idx.qmat] <- pars[idx.q]
  diag(qmat) <- -rowSums(qmat)
  qmat
}

starting.point.musse <- function(tree, k, q.div=5, yule=FALSE) {
  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  r <- if  ( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]
  p <- rep(c(pars.bd, r / q.div), c(k, k, k * (k-1)))
  names(p) <- argnames.musse(NULL, k)
  p
}
