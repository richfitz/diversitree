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
make.bd.nee <- function(tree, sampling.f=NULL, unresolved=NULL,
                        times=NULL) { 
  cache <- make.cache.bd.nee(tree, times, sampling.f, unresolved)

  ## This allows r < 0, a > 1.  a < 0 is not allowed though
  ## Also, if r < 0, then a must > 1 (and vv.).  XOR captures this.
  ##   if ( a < 0 || sign(1-a) != sign(r) ) return(-Inf)
  ## Alternatively, just enforce all(pars > 0)?
  ll <- function(pars, condition.surv=TRUE) {
    N <- cache$N
    x <- cache$x
    f <- cache$f
    unresolved <- cache$unresolved

    check.pars.bd(pars)
    if ( pars[2] == pars[1] )
      pars[1] <- pars[1] + 1e-12
    r <- pars[[1]] - pars[[2]]
    a <- pars[[2]] / pars[[1]]

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
  class(ll) <- c("bd", "function")
  ll
}

## 5: make.cache
make.cache.bd.nee <- function(tree=NULL, times=NULL,
                          sampling.f=NULL, unresolved=NULL) {
  if ( !is.null(times) && !is.null(tree) ) {
    stop("times cannot be specified if tree given")
  } else if ( is.null(times) && is.null(tree) ) {
    stop("Either times or tree must be specified")
  } else if ( is.null(times) ) {
    tree <- check.tree(tree)
    times <- as.numeric(sort(branching.times(tree), decreasing=TRUE))
  } else {
    times <- as.numeric(sort(times, decreasing=TRUE))
  }

  unresolved <- check.unresolved.bd(tree, unresolved)

  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, 1)

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
