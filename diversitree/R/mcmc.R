## Slice sampling

## These functions take a function 'f' that evaluates an argument 'x'
## which is a vector of parameters.  Across a single iteration,
## take the input x and return a vector of parameters 'y'
## corresponding to a new position.

## TODO: I want to allow arbitrary proposal methods.  Importantly, one
## useful feature will be to allow some integer positions to be
## updated with a slice sampler, and the others from a discrete
## distribution.  Hopefully this can be done via the 'proposal'
## vector.

## The default MCMC method will return a half-finished MCMC sample if
## interrupted with Ctrl-C.  It is expected that the default method
## will be sufficiently general to be useful for most approaches.
## Non-default methods should normally be default parameter setting.
## However, totally different MCMC algorithms could be used.

## Currently, the 'control' argument is not used by any method.
mcmc <- function(lik, x.init, nsteps, ...) {
  UseMethod("mcmc")
}

mcmc.default <- function(lik, x.init, nsteps, w, prior=NULL,
                         sampler=sampler.slice, fail.value=-Inf,
                         lower=-Inf, upper=Inf, print.every=1,
                         control=list(), ...) {
  npar <- length(x.init)
  if ( is.null(names(x.init)) )
    try(names(x.init) <- argnames(lik), silent=TRUE)
  
  if ( is.null(prior) )
    posterior <- protect(function(x) lik(x, ...),
                         fail.value.default=fail.value)
  else
    posterior <- protect(function(x) lik(x, ...) + prior(x, ...),
                         fail.value.default=fail.value)

  y.init <- posterior(x.init, fail.value=NULL)

  if ( !is.finite(y.init) || y.init == fail.value )
    stop("Starting point must have finite probability")

  check.par <- function(x) {
    if ( length(x) == 1 )
      rep(x, npar)
    else if ( length(x) == npar )
      x
    else
      stop(sprintf("'%s' of incorrect length",
                   deparse(substitute(x))))
  }

  lower <- check.par.length(lower, npar)
  upper <- check.par.length(upper, npar)
  w     <- check.par.length(w,     npar) # TODO: may need revisiting.

  check.bounds(lower, upper, x.init)

  hist <- vector("list", nsteps)

  if ( is.null(sampler) )
    sampler <- sampler.slice

  mcmc.loop <- function() {
    for ( i in seq_len(nsteps) ) {
      hist[[i]] <<- tmp <- sampler(posterior, x.init, y.init, w,
                                   lower, upper, control)
      x.init <- tmp[[1]]
      y.init <- tmp[[2]]
      if ( print.every > 0 && i %% print.every == 0 )
        cat(sprintf("%d: {%s} -> %2.5f\n", i,
                    paste(sprintf("%2.4f", tmp[[1]]), collapse=", "),
                    tmp[[2]]))
    }
    hist
  }

  mcmc.recover <- function(...) {
    hist <- hist[!sapply(hist, is.null)]
    if ( length(hist) == 0 )
      stop("MCMC was stopped before any samples taken")
    warning("MCMC was stopped prematurely: ", length(hist), "/", nsteps,
            " steps completed.  Truncated chain is being returned.",
            immediate.=TRUE)
    hist
  }

  hist <- tryCatch(mcmc.loop(), interrupt=mcmc.recover)
  hist <- cbind(i=seq_along(hist),
                as.data.frame(t(sapply(hist, unlist))))
  names(hist)[ncol(hist)] <- "p"

  hist
}

## This is common, so this helps reduce code duplication.
mcmc.lowerzero <- function(lik, x.init, nsteps, lower=0, ...)
  NextMethod("mcmc", lower=lower)

make.unipar <- function(f, x, i) {
  function(z) {
    x[i] <- z
    f(x)
  }
}
