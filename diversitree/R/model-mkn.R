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

make.mkn <- function(tree, states, k, strict=TRUE, control=list()) {
  control <- check.control.mkn(control, k)
  if ( control$method == "ode" )
    ll <- make.mkn.ode(tree, states, k, strict, control)
  else # exp or mk2
    ll <- make.mkn.exp(tree, states, k, strict, control)
  ll
}

make.mk2 <- function(tree, states, strict=TRUE, control=list()) {
  if ( !(is.null(control$method) || control$method == "mk2") )
    stop("Invalid control$method value (just omit it)")
  control$method <- "mk2"
  ll <- make.mkn(tree, states + 1, 2, strict, control)
  class(ll) <- c("mk2", "mkn", "function") # c("mk2", class(ll))
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

    base <- ceiling(log10(k + .5))
    fmt <- sprintf("q%%0%dd%%0%dd", base, base)
    sprintf(fmt, rep(1:k, each=k-1),
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
  if ( any(duplicated(value)) )
    stop("Duplicate argument names")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
## The use of nlminb is here for consistency with ape.  However, this
## is not a great optimisation method, and may cause problems for some
## cases.
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

## The other functions are all done by sub-models:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches

## Common utility functions for all Mkn-based models.
check.control.mkn <- function(control, k) {
  if ( is.null(control$method) ) {
    control$method <- "exp"
    if ( is.null(control$use.mk2) )
      control$use.mk2 <- FALSE
    else if ( k != 2 && control$use.mk2 )
      stop("control$use.mk2=TRUE only valid when k=2")
  } else if ( control$method == "mk2" ) {
    control$method <- "exp"
    control$use.mk2 <- TRUE
  } else if ( control$method == "exp" ) {
    if ( is.null(control$use.mk2) )
      control$use.mk2 <- FALSE
    else if ( k != 2 && control$use.mk2 )
      stop("control$use.mk2=TRUE only valid when k=2")
  } else if ( control$method != "ode" ) {
    stop("Unknown method")
  }
  control
}

check.pars.mkn <- function(pars, k) {
  if ( length(pars) != k*(k-1) )
    stop(sprintf("Invalid length parameters (expected %d)", k*(k-1)))
  if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
  TRUE
}

## Makes a function that converts k(k-1) parameters into a k^2 Q
## matrix.
make.mkn.pars <- function(k) {
  qmat <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))

  function(pars) {
    qmat[idx] <- pars
    diag(qmat) <- -rowSums(qmat)
    qmat
  }
}

## This is not the most efficient possible, but should not really be
## used in code for which this is a problem.  This just wraps around
## the above function so that this can be done as a once-off.
mkn.Q <- function(pars, k) {
  if ( missing(k) )
    k <- (1 + sqrt(1 + 4*(length(pars))))/2
  if ( abs(k - round(k)) > 1e-6 || length(pars) != k*(k-1) )
    stop("Invalid parameter length")
  make.mkn.pars(k)(pars)
}

## Additional functions:
stationary.freq.mkn <- function(pars) {
  if ( length(pars) == 2 )
    pars[2:1] / sum(pars)
  else
    .NotYetImplemented()
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
