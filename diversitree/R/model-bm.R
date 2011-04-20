## Brownian motion model from Luke's geiger package, so that I can use
## this without loading the entire package (and also have some more
## fun with different models in the future).

## This is actually slightly slower than Luke's version at the
## moment, which might be due to the shift from a log basis to a
## linear basis (the log basis is probably "easier", linearising the
## effects and helping optimisation).

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
## make.bm.old <- function(tree, states, meserr=NULL) {
##   cache <- make.cache.bm(tree, states, meserr)

##   n.tip <- cache$n.tip
##   states <- cache$states
##   meserr <- cache$meserr
##   vcv <- cache$vcv
##   one <- rep(1, n.tip)

##   if ( is.null(meserr) ) {
##     VI.tmp <- solve(vcv)
##     lik <- function(x) {
##       if ( x < 0 )
##         return(-Inf)
##       VI <- VI.tmp / x
##       ## By my calculations t(1) %*% VI %*% 1 = sum(VI)
##       ## t(one) %*% VI %*% states = sum(colSums(m) * x)
##       mu <- solve(t(one) %*% VI %*% one) %*% (t(one) %*% VI %*% states)
##       # from mvtnorm
##       dmvnorm(states, rep(mu, n.tip), vcv * x, log=TRUE)
##     }
##   } else {
##     lik <- function(x) {
##       if ( x < 0 )
##         return(-Inf)
##       vv <- x*vcv
##       if ( !is.null(meserr) )
##         diag(vv) <- diag(vv) + meserr^2
##       VI <- solve(vv)
##       mu <- solve(t(one) %*% VI %*% one) %*% t(one) %*% VI %*% states
##       dmvnorm(states, rep(mu, n.tip), vv, log=TRUE)
##     }
##   }
##   class(lik) <- c("bm", "function")
##   lik
## }

make.bm <- function(tree, states, meserr=NULL) {
  cache <- make.cache.bm(tree, states, meserr)

  n.tip <- cache$n.tip
  states <- cache$states
  meserr <- cache$meserr
  vcv <- cache$vcv
  one <- rep(1, n.tip)

  if ( is.null(meserr) ) {
    VI.tmp <- solve(vcv)
    lik <- function(x) {
      if ( x < 0 )
        return(-Inf)
      if ( length(x) != 1 )
        return(-Inf)
      VI <- VI.tmp / x
      ## By my calculations t(1) %*% VI %*% 1 = sum(VI)
      ## t(one) %*% VI %*% states = sum(colSums(VI) * states)
      ## mu <- solve(t(one) %*% VI %*% one) %*% (t(one) %*% VI %*% states)
      mu <- sum(colSums(VI) * states) / sum(VI)

      dmvnorm2(states, rep(mu, n.tip), vcv * x, VI, log=TRUE)
    }
  } else {
    lik <- function(x) {
      if ( x < 0 )
        return(-Inf)
      if ( length(x) != 1 )
        return(-Inf)
      vv <- x*vcv
      if ( !is.null(meserr) )
        diag(vv) <- diag(vv) + meserr^2
      VI <- solve(vv)
      mu <- sum(colSums(VI) * states) / sum(VI)
      dmvnorm2(states, rep(mu, n.tip), vv, VI, log=TRUE)
    }
  }
  class(lik) <- c("bm", "function")
  lik
}

## 2: print
print.bm <- function(x, ...) {
  cat("Brownian Motion likelihood function (1D)\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.bm <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    "beta" # Luke's name
  else
    ret
}
`argnames<-.bm` <- function(x, value) {
  if ( length(value) != 1 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x
}

## 4: find.mle
find.mle.bm <- function(func, x.init, method,
                         fail.value=NA, ...) {
  ## TODO: should use optimize, as this is a 1d problem.  It is
  ## bounded on [0, ?]; an evential wrapper could be written to run
  ## the bracketing.
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.bm")
}

mcmc.bm <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache
make.cache.bm <- function(tree, states, meserr) {
  tree <- check.tree(tree, ultrametric=FALSE)
  states <- check.states(tree, states)
  
  n.tip <- length(tree$tip.label)

  if ( is.null(meserr) )
    meserr <- NULL
  else if ( length(meserr) == 1 )
    meserr <- rep(meserr, n.tip)
  else
    meserr <- check.states(tree, meserr)
      
  list(tree=tree,
       states=states,
       meserr=meserr,
       n.tip=n.tip,
       vcv=vcv.phylo(tree))
}

dmvnorm2 <- function(x, mean, sigma, sigma.inv, log=FALSE) {
  distval <- mahalanobis(x, center=mean, cov=sigma.inv, TRUE)
  logdet <- as.numeric(determinant(sigma, TRUE)$modulus)
  logretval <- -(length(x) * log(2 * pi) + logdet + distval)/2
  if ( log )
    logretval
  else
    exp(logretval)
}
