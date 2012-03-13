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
make.bm <- function(tree, states, states.sd=0, control=list()) {
  control <- check.control.continuous(control)
  if ( control$method == "vcv" )
    make.bm.vcv(tree, states, states.sd)
  else
    make.bm.direct(tree, states, states.sd)
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
    "s2" # Luke's name was 'beta'
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
