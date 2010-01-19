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
  cache <- make.cache.bd(tree, times, sampling.f, unresolved)
  ll <- function(pars, ...) ll.bd(cache, c(pars, 0), ...)
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

  if ( !is.null(environment(func)$prior) )
    stop("Cannot yet get MAP point for Yule with prior - ",
         "use constrained bd model instead")

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
  else if ( is.null(times) )
    times <- as.numeric(sort(branching.times(tree), decreasing=TRUE))
  else
    times <- as.numeric(sort(times, decreasing=TRUE))
  
  if ( !is.null(sampling.f) ) .NotYetUsed("sampling.f")
  if ( !is.null(unresolved) ) .NotYetUsed("unresolved")

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

  list(N=N, x=x, tot.len=tot.len, n.node=n.node)
}

## make here only requires the ll function.  It cannot use cleanup
## though, as the normal cache structure is not generated.
ll.bd <- function(cache, pars, prior=NULL, condition.surv=TRUE) {
  N <- cache$N
  x <- cache$x

  r <- pars[1] - pars[2]
  a <- pars[2] / pars[1]
  ## This allows r < 0, a > 1.  a < 0 is not allowed though
  ## Also, if r < 0, then a must > 1 (and vv.).  XOR captures this.
  ##   if ( a < 0 || sign(1-a) != sign(r) ) return(-Inf)
  ## Alternatively, just enforce all(pars > 0)?
  if ( a < 0 || xor(a > 1, r < 0) )
    return(-Inf)

  if ( condition.surv )
    loglik <- lfactorial(N - 1) + (N - 2) * log(abs(r)) + r * sum(x[3:N]) +
      N * log(abs(1 - a)) - 2*sum(log(abs(exp(r * x[2:N]) - a)))
  else
    loglik <- lfactorial(N - 1) + (N - 2) * log(abs(r)) + r * sum(x[3:N]) +
      N * log(abs(1 - a)) - 2*sum(log(abs(exp(r * x[2:N]) - a))) +
        2*log(abs((1-a)/(1-a*exp(-r * x[2])))) + log(r/(1-a))


  ## Copied from cleanup
  if ( is.null(prior) )
    p <- loglik
  else if ( is.numeric(prior) )
    p <- loglik + prior.default(pars, prior)
  else if ( is.function(prior) )
    p <- loglik + prior(pars)
  else
    stop("Invalid 'prior' argument")

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
