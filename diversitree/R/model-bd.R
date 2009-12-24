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
make.bd <- function(tree, times=branching.times(tree),
                            sampling.f=NULL, unresolved=NULL) {
  if ( !is.null(sampling.f) ) .NotYetUsed("sampling.f")
  if ( !is.null(unresolved) ) .NotYetUsed("unresolved")

  if ( !missing(times) && !missing(tree) )
    stop("times cannot be specified if tree given")
  x <- c(NA, times)
  N <- length(x)
  cache <- list(N=N, x=x)

  ll <- function(pars, ...) bd.ll(cache, pars, ...)
  class(ll) <- c("bd", "function")
  ll
}

## 2: print
print.bd <- function(x, ...) {
  cat("Constant rate birth-death likelihood function:\n")
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

## 4: find.mle
find.mle.bd <- function(func, x.init, method,
                        fail.value=NA, ...) {
  ## I really should use parameters estimated from the Yule model
  ## here.  Currently using parameters from nowhere.
  if ( missing(x.init) )
    x.init <- structure(c(.2, .1), names=argnames(func))
  if ( missing(method) )
    method <- "nlm"
  NextMethod("find.mle", method=method, class.append="fit.mle.bd")
}

## make here only requires the ll function.  It cannot use cleanup
## though, as the normal cache structure is not generated.
bd.ll <- function(cache, pars, prior=NULL) {
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

  loglik <- lfactorial(N - 1) + (N - 2) * log(abs(r)) + r * sum(x[3:N]) +
    N * log(abs(1 - a)) - 2*sum(log(abs(exp(r * x[2:N]) - a)))

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
