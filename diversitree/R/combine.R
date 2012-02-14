## This comes from an email I wrote to Emma on 29 October 2010.
combine <- function(..., funs=list(...)) {
 if ( length(unique(lapply(funs, class))) != 1 )
   stop("All functions must have the same class")
 if ( length(unique(lapply(funs, argnames))) != 1 )
   stop("All functions must have the same argnames")
 
 ret <- function(pars, ...) {
   ans <- lapply(funs, function(f) f(pars, ...))
   sum(unlist(ans))
 }

 class(ret) <- c("combined", class(funs[[1]]))
 ret
}

argnames.combined <- function(x, ...)
  argnames(environment(x)$funs[[1]])

"argnames<-.combined" <- function(x, value) {
  argnames(environment(x)$funs[[1]]) <- value
  x
}

## This is even sweeter:
'+.bisse' <- function(x, y) combine(x, y)

## g <- liks[[1]]
## argnames(g) <- c("l0", "l1", "m0", "m1", "q0", "q1")

## f <- combine(funs=liks)
## p <- c(.1, .1, 0, 0, .03, .03)
## sum(sapply(liks, function(f) f(p)))

## argnames(f) <- c("l0", "l1", "m0", "m1", "q0", "q1")
## argnames(f)
## f(p)

## g <- constrain(f, l1 ~ l0)
## g(p[-2])


## It should then be possible to do this:

## chaparral.trees <- list(phy.cea, phy.arc)
## chaparral.states <- list(chars.cea, chars.arc)
## chaparral.samp <- list(c(0.913, 0.941, 0.875), c(0.674, 0.533, 0.75))

## lnL.each <- list()
## for (i in seq_along(chaparral.trees) )
##  lnL.each[[i]] <- make.geosse(chaparral.trees[[i]],
##                               chaparral.states[[i]],
##                               sampling.f = chaparral.samp[[i]])

## lik.multi <- constrain(combine(lnL.each), sAB ~ 0)

## You could possibly also overload '+' for likelihood functions (argument names incorrect, I think)
## '+.geosse' <- function(x, y) combine(list(x, y))

## And do it this way
##  lik1 <- make.geosse(...whatever...)
##  lik2 <- make.geosse(...whatever...)
##  lik.multi <- lik1 + lik2

