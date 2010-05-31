## The B-D model is different to the others in that I am not using
## most of the infrastructure - instead the entire calculations are
## done at once.

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
make.bd <- function(tree, times=NULL,
                    sampling.f=NULL, unresolved=NULL) {
  cache <- make.cache.bd(tree, times, sampling.f, unresolved)
  ll <- function(pars, ...) ll.bd(cache, pars, ...)
  class(ll) <- c("bd", "function")
  ll
}

make.yule <- function(tree, times=NULL,
                      sampling.f=NULL, unresolved=NULL) {
  if ( !is.null(sampling.f) || !is.null(unresolved) )
    stop("Cannot yet do Yule model with unresolved clades or sampling.f")
  cache <- make.cache.bd(tree, times, sampling.f, unresolved)
  ll <- function(pars, ...)
    ll.bd(cache, c(pars, 0), ...)
  class(ll) <- c("yule", "bd", "function")
  ll
}

## 2: print
print.bd <- function(x, ...) {
  cat("Constant rate birth-death likelihood function:\n")
  print(unclass(x))
}

print.yule <- function(x, ...) {
  cat("Constant rate Yule model likelihood function:\n")
  print(unclass(x))
}


## 3: argnames / argnames<-
argnames.bd <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("lambda", "mu")
  else
    ret
}
`argnames<-.bd` <- function(x, value) {
  if ( length(value) != 2 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

argnames.yule <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    "lambda"
  else
    ret
}
`argnames<-.yule` <- function(x, value) {
  if ( length(value) != 1 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.bd <- function(func, x.init, method,
                        fail.value=NA, ...) {
  ## I really should use parameters estimated from the Yule model
  ## here.  Currently using parameters from nowhere, as doing this is
  ## slightly less trivial than one would hope (need to get access to
  ## the cache, which might be hidden by a constraint...
  if ( missing(x.init) ) {
    ## tmp <- find.mle.yule(func, ...)
    ## x.init <- structure(c(tmp$par, 0), names=argnames(func))
    warning("Guessing initial parameters - may do badly")
    x.init <- structure(c(.2, .1), names=argnames(func))
  }
  if ( missing(method) )
    method <- "nlm"
  find.mle.default(func, x.init, method, fail.value,
                   class.append="fit.mle.bd", ...)
  ## NextMethod("find.mle", x.init=x.init, method=method,
  ##           class.append="fit.mle.bd")
}

find.mle.yule <- function(func, x.init, method, fail.value=NA,
                          ...) {
  ## TODO: I think this will fail after constraining, which should not
  ## be done.
  ## TODO: This is only true for certain cases (no funny business with
  ## sampling.f/unresolved - this might require moving to
  ## numeric solutions?)
  cache <- environment(func)$cache
  condition.surv <- list(...)$condition.surv
  if ( is.null(condition.surv) )
    condition.surv <- TRUE

  n.node <- if ( condition.surv ) cache$n.node - 1 else cache$n.node
  lambda <- n.node / cache$tot.len
  obj <- list(par=c(lambda=lambda),
              lnLik=func(lambda, ...),
              counts=NA,
              code=0,
              gradient=NA,
              method="analytic")
  class(obj) <- c("fit.mle.yule", "fit.mle")
  obj
}

## 5: make.cache
make.cache.bd <- function(tree=NULL, times=NULL,
                          sampling.f=NULL, unresolved=NULL) {
  if ( !is.null(times) && !is.null(tree) )
    stop("times cannot be specified if tree given")
  else if ( is.null(times) && is.null(tree) )
    stop("Either times or tree must be specified")
  else if ( is.null(times) ) {
    tree <- check.tree(tree)
    times <- as.numeric(sort(branching.times(tree), decreasing=TRUE))
  } else
    times <- as.numeric(sort(times, decreasing=TRUE))

  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, 1)

  if ( !is.null(unresolved) && length(unresolved) > 0 ) {
    if ( is.null(names(unresolved)) || !is.numeric(unresolved) )
      stop("'unresolved' must be a named numeric vector")
    if ( !(all(names(unresolved) %in% tree$tip.label)) )
      stop("Unknown species in 'unresolved'")
    if ( any(unresolved < 0) )
      stop("All unresolved entries must be > 0")

    if ( any(unresolved == 0) ) {
      ## TODO: surely this should be entries that are one?
      warning("Removing unresolved entries that are zero")
      unresolved <- unresolved[unresolved != 0]
    }

    if ( length(unresolved) == 0 )
      unresolved <- NULL
    else {
      i <- match(names(unresolved), tree$tip.label)
      unresolved <- list(n=unresolved,
                         t=tree$edge.length[match(i, tree$edge[,2])])
    }
  } else {
    unresolved <- NULL
  }

  x <- c(NA, times)
  N <- length(x)

  if ( is.null(tree) ) {
    tot.len <- sum(diff(times[1] - c(times, 0)) *
                   2:(length(times) + 1))
    n.node <- length(times)
  } else {
    tot.len <- sum(tree$edge.length)
    n.node <- tree$Nnode
  }

  list(N=N, x=x, f=sampling.f, unresolved=unresolved, tot.len=tot.len,
       n.node=n.node)
}

## This allows r < 0, a > 1.  a < 0 is not allowed though
## Also, if r < 0, then a must > 1 (and vv.).  XOR captures this.
##   if ( a < 0 || sign(1-a) != sign(r) ) return(-Inf)
## Alternatively, just enforce all(pars > 0)?
ll.bd <- function(cache, pars, prior=NULL, condition.surv=TRUE) {
  if ( !is.null(prior) )
    stop("'prior' argument to likelihood function no longer accepted")

  N <- cache$N
  x <- cache$x
  f <- cache$f
  unresolved <- cache$unresolved

  if ( length(pars) != 2 )
    stop("Incorrect number of parameters")
  if ( pars[2] == pars[1] )
    pars[1] <- pars[1] + 1e-12

  r <- pars[1] - pars[2]
  a <- pars[2] / pars[1]

  if ( pars[1] <= 0 || pars[2] < 0 )
    return(-Inf)

  if ( f < 1 )
    loglik <- lfactorial(N - 1) +
      (N - 2) * log(f*abs(r)) + N * log(abs(1 - a)) +
        r * sum(x[3:N]) -
          2*sum(log(abs(f * exp(r * x[2:N]) - a + 1 - f)))
  else # this may not be faster enough to warrant keeping.
    loglik <- lfactorial(N - 1) +
      (N - 2) * log(abs(r)) + N * log(abs(1 - a)) + r * sum(x[3:N]) -
        2*sum(log(abs(exp(r * x[2:N]) - a)))

  if ( !condition.surv )
    loglik <- loglik +
      log(f * f * r * (1 - a)) -
        2*log(abs(exp(-r * x[2])*(a - 1 + f) - f))

  if ( !is.null(unresolved) ) {
    ert <- exp(r * unresolved$t)
    loglik <- loglik +
      sum((unresolved$n-1) * (log(abs(ert - 1)) - log(abs(ert - a))))
  }

  loglik
}

starting.point.bd <- function(tree, yule=FALSE) {
  pars.yule <- c(coef(find.mle(make.yule(tree))), 0)
  if ( yule )
    p <- pars.yule
  else
    p <- coef(find.mle(make.bd(tree), pars.yule))
  names(p) <- c("lambda", "mu")
  p
}
