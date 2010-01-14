## The Mk2 model of character evolution.  This is likely to end up
## being the Mkn model, really.  However, a simple 'mk2' interface
## will remain that does not require specifying the entire transition
## matrix.

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
##   9. branches.unresolved

## 1: make
make.mk2 <- function(tree, states) {
  cache <- make.cache.mkn(tree, states + 1, 2)
  ll <- function(pars, ...) ll.mk2(cache, pars, ...)
  class(ll) <- c("mk2", "mkn", "function")
  attr(ll, "k") <- 2
  ll
}

make.mkn <- function(tree, states, k) {
  cache <- make.cache.mkn(tree, states, k)
  branches <- make.branches.mkn(k)
  qmat <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))

  ll.mkn <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     root.p=NULL, intermediates=FALSE) {
    if ( length(pars) != k*(k-1) )
      stop("Invalid length parameters")
    if ( any(!is.finite(pars)) )
      return(-Inf)
    qmat[idx] <- pars
    diag(qmat) <- -rowSums(qmat)
    
    ans <- all.branches(qmat, cache, initial.conditions.mkn, branches,
                        branches.unresolved.mkn)
    loglik <- root.mkn(ans$init[cache$root,], ans$lq, pars, root, root.p)
    cleanup(loglik, pars, prior, intermediates, cache, ans)
    ##     if ( !is.null(prior) )
    ##       loglik <- loglik + prior.mkn(pars, prior)
    ##     loglik
  }

  ll <- function(pars, ...) ll.mkn(cache, pars, ...)
  class(ll) <- c("mkn", "function")
  attr(ll, "k") <- k  
  ll
}

## 2: print
print.mkn <- function(x, ...) {
  cat("Mk-n likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.mkn <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    k <- attr(x, "k")
    sprintf("q%d%d", rep(1:k, each=k-1),
            unlist(lapply(1:k, function(i) (1:k)[-i])))
  } else {
    ret
  }
}

argnames.mkn <- function(x, k, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    if ( missing(k) )
      k <- attr(x, "k")
    else
      if ( !is.null(x) )
        stop("k can only be be given if x is null")
    sprintf("q%d%d", rep(1:k, each=k-1),
            unlist(lapply(1:k, function(i) (1:k)[-i])))
  } else {
    ret
  }
}
argnames.mk2 <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("q01", "q10")
  else
    ret
}
`argnames<-.mkn` <- function(x, value) {
  k <- environment(x)$cache$k  
  if ( length(value) != k*(k-1) )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.mkn <- function(func, x.init, method,
                         fail.value=NA, ...) {
  ## These parameters are just pulled out of thin air
  if ( missing(x.init) )
    x.init <- structure(c(.1, .1), names=argnames(func))
  if ( missing(method) )
    method <- "nlminb"
  NextMethod("find.mle", method=method, class.append="fit.mle.mkn")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.mkn <- function(tree, states, k=2) {
  if ( is.null(names(states)) )
    stop("The states vector must contain names")
  if ( !all(tree$tip.label %in% names(states)) )
    stop("Not all species have state information")
  states <- states[tree$tip.label]
  names(states) <- tree$tip.label
  cache <- make.cache(tree)
  cache$k <- k
  cache$tip.state  <- states
  cache$y <- initial.tip.mkn(cache)
  cache
}
initial.tip.mkn <- function(cache) {
  k <- cache$k
  tip.state <- cache$tip.state
  if ( any(tip.state < 1 | tip.state > k, na.rm=TRUE) )
    stop(sprintf("tip states must be in the range [1, %d]", k))
  
  y <- matrix(rep(rep(0, k), k + 1), k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,1:k]) <- 1
  i <- cache$tip.state
  i[is.na(i)] <- k + 1
  list(y=y, i=i, types=sort(unique(i)))
}
root.mkn <- function(vals, lq, pars, root, root.p=NULL) {
  k <- length(vals)
  if ( !is.null(root.p) ) {
    if ( root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    else if ( length(root.p) != k )
      stop("root.p of wrong length")
  }

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    p <- stationary.freq.mkn(pars)
  else if ( root == ROOT.OBS )
    p <- vals/sum(vals)
  else if ( root == ROOT.GIVEN )
    p <- root.p
  else if ( root != ROOT.BOTH )
    stop("Invalid root mode")

  logcomp <- sum(lq)
  if ( root == ROOT.BOTH )
    loglik <- log(vals) + logcomp
  else
    loglik <- log(sum(p * vals)) + logcomp
  loglik
}

## 6: ll (ll.mkn is done within make.mkn)
ll.mk2 <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                   root.p=NULL, intermediates=FALSE) {
  if ( length(pars) != 2 )
    stop("Invalid length parameters")
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)
  ans <- all.branches(pars, cache, initial.conditions.mkn, branches.mk2,
                      branches.unresolved.mkn)
  loglik <- root.mkn(ans$init[cache$root,], ans$lq, pars, root, root.p)
  ##   if ( !is.null(prior) )
  ##     loglik <- loglik + prior.mkn(pars, prior)
  ##   loglik
  cleanup(loglik, pars, prior, intermediates, cache, ans)
}

## 7: initial.conditions:
initial.conditions.mkn <- function(init, pars, t, is.root=FALSE)
  init[1,] * init[2,]

## 8: branches (separate for mk2 and mkn)
branches.mk2 <- function(y, len, pars, t0) {
  q01 <- pars[1]
  q10 <- pars[2]
  if ( q01 + q10 > 0 ) {
    x <- exp(-(q01+q10)*len) * (y[1] - y[2])
    z <- q10 * y[1] + q01 * y[2]
    ret <- cbind(z + x * q01, z - x * q10) / (q01 + q10)
  } else { # Special case...
    ret <- matrix(rep(y, length(len)), length(len), 2, TRUE)
  }

  q <- rowSums(ret)
  i <- q > 0
  ret[i,] <- ret[i,] / q[i]
  lq <- q
  lq[i] <- log(q[i])
  cbind(lq, ret, deparse.level=0)
}

## The n-state version is not much different:
make.branches.mkn <- function(k) {
  RTOL <- ATOL <- 1e-8
  eps <- 0
  
  if ( k == 2 )
    warning("Two states is faster with Mk2")
  mkn.ode <- make.ode("derivs_mkn", "diversitree", "initmod_mkn", k, FALSE)

  branches.mkn <- function(y, len, pars, t0)
    t(mkn.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  make.branches(branches.mkn, 1:k)
  
##   function(y, len, pars, t0) {
##     ret <- t(mkn.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
##     if ( all(ret >= eps) ) {
##       q <- apply(ret, 1, min)
##       i <- q > 0
##       ret[i,] <- ret[i,] / q[i]
##       lq <- q
##       lq[i] <- log(q[i])
##       cbind(lq, ret, deparse.level=0)
##     } else {
##       ti <- len[length(len)]/2
##       len1 <- c(len[len <= ti], ti)
##       len2 <- len[len > ti] - ti
##       n1 <- length(len1)
##       ret1 <- Recall(y, len1, pars, t0)
##       ret2 <- Recall(ret1[n1,-1], len2, pars, t0 + ti)
##       ret2[,1] <- ret2[,1] + ret1[n1,1]
##       rbind(ret1[-n1,], ret2)
##     }
##   }
}

## 9: branches.unresolved
branches.unresolved.mkn <- function(...)
  stop("Cannot use unresolved clades with Mk2/Mkn")

## Additional functions:
stationary.freq.mkn <- function(pars) {
  .NotYetImplemented()
}


