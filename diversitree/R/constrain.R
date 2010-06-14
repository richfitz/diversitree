## This is still useful:
argnames <- function(x, ...)
  UseMethod("argnames")
`argnames<-` <- function(x, value)
  UseMethod("argnames<-")
## argnames.default <- function(x, ...)
##   attr(x, "argnames")
## `argnames<-.default` <- function(x, value) {
##   attr(x, "argnames") <- value
##   x
## }
`argnames<-.constrained` <- function(x, value) {
  stop("Cannot set argnames on constrained function")
}
argnames.constrained <- function(x, ...)
  environment(x)$names.rhs


## The LHS of a formula must be a single variable name that exists in
## "names.lhs"
##
## The RHS can be one of
##   - numeric value
##   - expression
## If it is an expression, then all variable names must be found in
## names.rhs (or perhaps in the containing environment - check in the
## future?)

## I might eventually allow formulae of the form
##   lambda | lambda1 ~ lambda0
## to allow renaming?

## How do I use this to allow recasting to alternative bases?
## Forward definitions for just the diversification rate are
##   foo(f, r0 ~ lambda0 - mu0, r1 ~ lambda1 - mu1)
## and for both diversification and relative extinction    
##   foo(f,
##       r0 ~ lambda0 - mu0, r1 ~ lambda1 - mu1,
##       e0 ~ mu0 / lambda0, e1 ~ mu1 / lambda1)
## Backward definitions (leaving mu unchanged)
##   foo(f, lambda0 ~ r0 + mu0, lambda1 ~ r1 + mu1)
## both:
##   foo(f, lambda0 ~ r0/(1 - e0), lambda1 ~ r1/(1 - e1),
##       r0 * e0 / (1 - e0), r1 * e1 / (1 - e1))

## TODO: Check that lhs does not appear on the rhs
constrain.parse <- function(formula, names.lhs, names.rhs,
                            extra=NULL) {
  formula <- as.formula(formula)
  if ( length(formula) != 3L )
    stop("Invalid formula")
  lhs <- formula[[2]]
  rhs <- formula[[3]]

  ## Checking the lhs is easy: is the lhs in the list of allowable
  ## names and of length 1?  Anything that does not match this is
  ## invalid.
  if ( !is.name(lhs) || is.na(match(as.character(lhs), names.lhs)) )
    stop("Invalid target on LHS of formula" )

  ## Checking the rhs is more difficult.  We are OK if any of the
  ## following is met:
  ##   Numeric values (checked at the end)
  ##   If all vars in rhs are in names.rhs
  ##   There is a single var and it is in extra
  ## Failing that, if the rhs is a single variable that does exist in
  ## the calling environment.
  if ( is.language(rhs) ) {
    vars <- all.vars(rhs)
    ok <- (all(vars %in% names.rhs) ||
           length(vars) == 1 && vars %in% extra)
    if ( !ok && length(vars) == 1 ) {
      e <- parent.frame()
      if ( exists(vars, e) ) {
        rhs <- get(vars, e)
        ok <- TRUE
      }
    }

    if ( !ok )
      stop("Invalid RHS of formula:\n\t", as.character(rhs))
    if ( as.character(lhs) %in% vars )
      stop("LHS cannot appear in RHS")
  } else if ( !is.numeric(rhs) ) {
    stop("RHS must be expression, variable or number")
  }
  list(lhs, rhs)
}

## First up, consider the one-shot case: don't worry about incremental
## updates.

## For the first case, everything is OK on the lhs and rhs
## For subsequent cases:
##   lhs cannot contain things that are
##      - constrained things (already lhs anywhere)
##      - things constrained to (things on the rhs anywhere)
##   rhs cannot contain things that are
##      - constrained things (already lhs anywhere)
## It is possibly worth pulling out all the numerical constants and
## the "paired" parameters here to avoid using eval where
## unnecessary.  However, this makes the function substantially uglier
## for a very minor speedup.

constrain <- function(f, ..., formulae=NULL, names=argnames(f),
                      extra=NULL) {
  if ( inherits(f, "constrained") ) {
    formulae <- c(attr(f, "formulae"), formulae)
    f <- attr(f, "func")
  }

  formulae <- c(formulae, list(...))
  names.lhs <- names.rhs <- names
  rels <- list()
  
  for ( formula in formulae ) {
    res <- constrain.parse(formula, names.lhs, names.rhs, extra)
    names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
    names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
    rels <- c(rels, structure(res[2], names=as.character(res[[1]])))
  }

  i <- match(unique(sapply(rels, as.character)), extra)
  final <- c(extra[sort(i[!is.na(i)])], names.rhs)
  npar <- length(final)

  ## "free" are the parameters that have nothing special on their RHS
  ## and are therefore passed directly through
  free <- setdiff(names.rhs, names(rels))
  free.i <- match(free, names) # index in full variables
  free.j <- match(free, final) # index in given variables.

  ## Targets are processed in the same order as given by formulae. 
  target.i <- match(names(rels), names)

  pars <- rep(NA, length(names))
  names(pars) <- names
  g <- function(x, ..., pars.only=FALSE) {
    if ( length(x) != npar )
      stop(sprintf("x of wrong length (expected %d)", npar))

    pars[free.i] <- x[free.j]
    e <- structure(as.list(x), names=final)
    pars[target.i] <- unlist(lapply(rels, eval, e))

    if ( pars.only )
      pars
    else
      f(pars, ...)
  }

  class(g) <- c("constrained", class(f))
  attr(g, "argnames") <- final
  attr(g, "formulae") <- formulae
  attr(g, "extra") <- extra
  attr(g, "func") <- f
  g
}

argnames.constrained <- function(x, ...)
  attr(x, "argnames")

print.constrained <- function(x, ...) {
  rels <- environment(x)$rels
  NextMethod("print")
  cat("Free parameters:", paste(argnames(x), collapse=", "), "\n")
  cat("Constraints:\n")
  cat(sprintf("\t%s ~ %s\n", names(rels),
              unlist(lapply(rels, deparse))))
}

update.constrained <- function(object, free, ...) {
  func <- attr(object, "func")
  formulae <- attr(object, "formulae")
  extra <- attr(object, "extra")

  lhs <- sapply(formulae, function(x) as.character(x[[2]]))
  i <- match(free, lhs)
  if ( any(is.na(i)) )
    stop("Variable(s) ", paste(free[i], collapse=", "),
         " are not constrained")
  formulae <- formulae[-i]
  if ( length(formulae) > 0 )
    constrain(func, formulae=formulae, extra=extra)
  else
    func
}
