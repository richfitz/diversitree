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

## This version does not use all.branches, but is substantially faster
## (c.f. model-mkn-legacy, left for didactic purposes).
make.mk2 <- function(tree, states) {
  ll <- make.mkn(tree, states + 1, 2, TRUE)
  class(ll) <- c("mk2", "mkn", "function")
  ll
}
make.mkn <- function(tree, states, k, use.mk2=FALSE) {
  cache <- make.cache.mkn(tree, states, k, use.mk2)

  qmat <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))

  len.uniq <- cache$len.uniq
  len.idx <- cache$len.idx
  n.tip <- cache$n.tip

  ll.mkn <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     root.p=NULL, intermediates=FALSE) {
    if ( !is.null(prior) )
      stop("'prior' argument to likelihood function no longer accepted")
    if ( length(pars) != k*(k-1) )
      stop(sprintf("Invalid length parameters (expected %d)", k*(k-1)))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    qmat[idx] <- pars
    diag(qmat) <- -rowSums(qmat)
    ans <- all.branches.mkn(qmat, cache)
    d.root <- ans$init[[cache$root]]
    root.p <- root.p.mkn(d.root, pars, root, root.p)
    loglik <- root.mkn(d.root, ans$lq, root.p)
    if ( intermediates ) {
      ans$init[seq_len(n.tip)] <- matrix.to.list(cache$y$y[cache$y$i,])
      ans$root.p <- root.p
    }

    cleanup(loglik, pars, intermediates, cache, ans)
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

mcmc.mkn <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.mkn <- function(tree, states, k, use.mk2) {
  tree <- check.tree(tree)
  states <- check.states(tree, states)

  cache <- make.cache(tree)
  cache$k <- k
  cache$tip.state  <- states
  ## TODO: deal with unknown states; tip calculations for unknown tip
  ## states just return the sum of nonzero ys from the matrix
  ## multiplication, so they can just be done separately.
  cache$map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
  cache$idx.tip <- cbind(c(cache$map[cache$tip.state,]),
                         rep(seq_len(cache$n.tip), k))
  cache$len.uniq <- sort(unique(cache$len))
  cache$len.idx <- match(cache$len, cache$len.uniq)
  cache$y <- initial.tip.mkn(cache) 

  ## Alter things to make it more speedy.  The '0' denotes base zero
  ## indices.
  cache$children0 <- as.integer(t(cache$children-1))
  cache$order0 <- as.integer(cache$order-1)
  cache$use.mk2 <- use.mk2
  cache$f.pij <- if ( k == 2 && use.mk2 )
    pij.mk2 else make.pij.mkn(k)
  
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

root.p.mkn <- function(vals, pars, root, root.p=NULL) {
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
  else if ( root == ROOT.BOTH )
    p <- NULL
  else
    stop("Invalid root mode")

  p
}

root.mkn <- function(vals, lq, root.p) {
  logcomp <- sum(lq)
  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(vals) + logcomp
  else
    loglik <- log(sum(root.p * vals)) + logcomp
  loglik
}

## 6: ll (done within make.mkn)

## 7: initial.conditions:
initial.conditions.mkn <- function(init, pars, t, is.root=FALSE)
  init[[1]] * init[[2]]

## 8: branches (separate for mk2 and mkn)
pij.mk2 <- function(len, pars) {
  ## The 3,2 indices here are because pars is a Q matrix by the time
  ## this gets called.
  q01 <- pars[3]
  q10 <- pars[2]
  x <- exp(-(q01+q10)*len)
  rbind((x*q01 + q10),
        (1 - x)*q10,
        (1 - x)*q01,
        (x*q10 + q01)) / (q01 + q10)
}

make.pij.mkn <- function(k) {
  pij.ode <- make.ode("derivs_mkn_pij", "diversitree",
                      "initmod_mkn", k*k, FALSE)
  ATOL <- RTOL <- 1e-8
  yi <- diag(k)
  
  function(len, pars)
    pij.ode(yi, c(0, len), pars,
            rtol=RTOL, atol=ATOL)[-1,-1,drop=FALSE]
}

## 9: branches.unresolved: n/a, as no option for it

## Additional functions:
stationary.freq.mkn <- function(pars) {
  if ( length(pars) == 2 )
    pars[2:1] / sum(pars)
  else
    .NotYetImplemented()
}
mkn.Q <- function(pars, k) {
  if ( missing(k) )
    k <- (1 + sqrt(1 + 4*(length(pars))))/2
  if ( abs(k - round(k)) > 1e-6 || length(pars) != k*(k-1) )
    stop("Invalid parameter length")
  qmat <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  qmat[idx] <- pars
  diag(qmat) <- -rowSums(qmat)
  qmat
}

## The new integrator.
all.branches.mkn <- function(pars, cache) {
  ## At this point, the parameters are assumed to be a Q matrix
  k <- cache$k

  idx.tip <- cache$idx.tip
  n.tip <- cache$n.tip
  n <- length(cache$len)
  order0 <- cache$order0

  len.uniq <- cache$len.uniq
  len.idx <- cache$len.idx
  f.pij <- cache$f.pij
  pij <- f.pij(len.uniq, pars)[,len.idx]
  
  lq <- numeric(n)
  branch.init <- branch.base <- matrix(NA, k, n)
  storage.mode(branch.init) <- "numeric"

  ## tips
  ans <- matrix(pij[idx.tip], n.tip, k)
  q <- rowSums(ans)
  branch.base[,seq_len(n.tip)] <- t(ans/q)
  lq[seq_len(n.tip)] <- log(q)

  ## branches
  ans <- .C("r_mkn_core",
            k        = as.integer(k),
            n        = length(order0) - 1L,
            order    = order0,
            children = cache$children0,
            pij      = pij,
            init     = branch.init,
            base     = branch.base,
            lq       = lq,
            NAOK=TRUE, DUP=FALSE)

  list(init=matrix.to.list(t(ans$init)),
       base=matrix.to.list(t(ans$base)),
       lq=ans$lq, pij=pij)
}

