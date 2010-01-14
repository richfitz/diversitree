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
constrain.parse <- function(formula, names.lhs, names.rhs) {
  formula <- as.formula(formula)
  if ( length(formula) != 3L )
    stop("Invalid formula")
  lhs <- formula[[2]]
  rhs <- formula[[3]]
  
  if ( !is.name(lhs) || is.na(match(as.character(lhs), names.lhs)) )
    stop("Invalid target on LHS of formula" )
  if ( is.language(rhs) ) {
    vars <- all.vars(rhs)
    if ( !all(vars %in% names.rhs) ) {
      if ( length(vars) == 1 && exists(vars) )
        ## TODO: Check that 'vars' is not really already in the
        ## function names.
        rhs <- get(vars)
      else
        stop("Some elements of the RHS were not found in names.rhs:\n\t",
             paste(setdiff(vars, names.rhs), collapse=", "))
    }
    if ( as.character(lhs) %in% vars )
      stop("LHS cannot appear in RHS")
  } else if ( !is.numeric(rhs) ) {
    stop("RHS must be expression, variable or number")
  }
  list(lhs, rhs)
}

## First up, consider the one-shot case: don't worry about incremental
## updates.

## OK, that's awesome.
## constrain.parse(lambda0 ~ lambda1, "lambda0", "lambda1")
## constrain.parse(lambda0 ~ lambda1 + mu0, "lambda0", "lambda1")

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
constrain <- function(f, ..., formulae=NULL, names=argnames(f)) {
  if ( inherits(f, "constrained") )
    warning("It is probaably not a good idea to constrain a constrained function")
  formulae <- c(list(...), formulae)
  names.lhs <- names.rhs <- names
  rels <- list()
  for ( formula in formulae ) {
    res <- constrain.parse(formula, names.lhs, names.rhs)
    names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
    names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
    rels <- c(rels, structure(res[2], names=as.character(res[[1]])))
  }

  npar <- length(names.rhs)
  free.i <- match(names.rhs, names)
  target.i <- match(names(rels), names)
  pars <- rep(NA, length(names))
  names(pars) <- names

  g <- function(x, ..., pars.only=FALSE) {
    if ( length(x) != npar )
      stop("x of wrong length")
    pars[free.i] <- x
    pars[target.i] <- unlist(lapply(rels, eval,as.list(pars[free.i])))
    if ( pars.only )
      pars
    else
      f(pars, ...)
  }
  class(g) <- c("constrained", class(f))
  g
}

print.constrained <- function(x, ...) {
  rels <- environment(x)$rels
  NextMethod("print")
  cat("Free parameters:", paste(argnames(x), collapse=", "), "\n")
  cat("Constraints:\n")
  cat(sprintf("\t%s ~ %s\n", names(rels),
              unlist(lapply(rels, deparse))))
}

