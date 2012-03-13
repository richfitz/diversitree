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

## At some point there will be 
make.bd <- function(tree, sampling.f=NULL, unresolved=NULL,
                    times=NULL, control=list()) {
  control <- check.control.bd(control, times)
  if ( control$method == "nee" )
    make.bd.nee(tree, sampling.f, unresolved, times)
  else
    make.bd.ode(tree, sampling.f, unresolved, control)
}

## Yule is somewhat weird, as we allow likelihood calculations, but
## cheat on the ML search and go straight for the ML point.
make.yule <- function(tree, sampling.f=NULL, unresolved=NULL,
                      times=NULL, control=list()) {
  control <- check.control.bd(control)
  if ( !is.null(sampling.f) || !is.null(unresolved) )
    stop("Cannot yet do Yule model with unresolved clades or sampling.f")
  ll.bd <- make.bd(tree, sampling.f, unresolved, times, control)
  ll <- function(pars, condition.surv=TRUE)
    ll.bd(c(pars, 0), condition.surv)
  
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
    warning("Guessing initial parameters - may do badly")
    x.init <- structure(c(.2, .1), names=argnames(func))
  }
  if ( missing(method) )
    method <- if (inherits(func, "bd.split")) "subplex" else "nlm"
  find.mle.default(func, x.init, method, fail.value,
                   class.append="fit.mle.bd", ...)
}

find.mle.yule <- function(func, x.init, method, fail.value=NA,
                          ...) {
  cache <- environment(environment(func)$ll.bd)$cache

  ## This analytic solution is only correct when the tree is fully
  ## resolved.  This is enforced by the make.yule function.
  condition.surv <- list(...)$condition.surv
  if ( is.null(condition.surv) )
    condition.surv <- TRUE

  n.node <- if ( condition.surv ) cache$n.node - 1 else cache$n.node
  lambda <- n.node / cache$tot.len
  obj <- list(par=c(lambda=lambda),
              lnLik=func(lambda, condition.surv),
              counts=NA,
              code=0,
              gradient=NA,
              method="analytic")

  ## This class here is needed so that this can be compared against a
  ## BD fit.
  class(obj) <- c("fit.mle.bd", "fit.mle")
  obj
}

mcmc.bd <- mcmc.lowerzero

starting.point.bd <- function(tree, yule=FALSE) {
  pars.yule <- c(coef(find.mle(make.yule(tree))), 0)
  if ( yule )
    p <- pars.yule
  else
    p <- coef(find.mle(make.bd(tree), pars.yule, method="subplex"))
  names(p) <- c("lambda", "mu")
  p
}

check.control.bd <- function(control, times) {
  control <- modifyList(list(method="nee"), control)
  if ( !(control$method %in% c("nee", "ode")) )
    stop(sprintf("Unknown method '%s'", control$method))
  control
}

check.pars.bd <- function(pars) {
  if ( length(pars) != 2 )
    stop("Invalid parameter length (expected 2)")
  if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
}
