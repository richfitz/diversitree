make.ode.gslode <- function(info, control) {
  check.gslode(error=TRUE)
 
  n.var <- info$ny
  n.par <- info$np 
  rtol <- atol <- control$tol
  stepper <- control$gsl.stepper

  ## Some checking:
  if ( length(rtol) != 1 )
    stop("rtol must (currently) be scalar")
  if ( length(atol) != 1 )
    stop("atol must (currently) be scalar")

  if ( is.function(info$derivs) ) {
    derivs <- info$derivs
    ode <- new(GslOdeR, derivs, environment(derivs), n.var)    
  } else {
    model <- info$name.ode
    dll <- info$dll
    derivs <- sprintf("derivs_%s_gslode", model)
    derivs <- getNativeSymbolInfo(derivs, PACKAGE=dll)$address

    if ( isTRUE(info$time.varying) ) {
      ode <- new(GslOdeTime, derivs, n.var, info$tm)
    } else {
      ode <- new(GslOdeCompiled, derivs, n.var)
    }
  }

  ## Control parameters (will get tweaked.
  ode$set_control(list(atol=atol, rtol=rtol, algorithm=stepper))
  
  function(vars, times, pars) {
    if ( length(pars) != n.par )
      stop("Incorrect parameter length")
    if ( length(vars) != n.var )
      stop("Incorrect variable length")
    if ( length(times) <= 1 )
      stop("Need >= 2 times")

    ## TODO: Will need to think about whether to drop first row
    ## always.
    ode$run(times, vars, pars)
  }
}

## This might not need checking, actually...
check.gslode <- function(error=TRUE) {
  ok <- exists("GslOdeCompiled")
  if ( error && !ok )
    stop("diversitree built without GslOde support")
  ok
}
