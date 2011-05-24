## This is a *very* basic interface to the the CVODES functions.

## TODO: Some documentation.

## Create a CVODES function.  This is similar in spirit to make.ode(),
## though the function definition of the returned functions differs.
cvodes <- function(n.var, n.par, derivs, rtol, atol, package=NULL) {
  check.cvodes(error=TRUE)
 
  n.var <- as.integer(n.var)
  n.par <- as.integer(n.par)
  if ( is.null(package) ) package <- ""
  derivs <- getNativeSymbolInfo(derivs, PACKAGE=package)$address
  rtol <- as.numeric(rtol)
  atol <- as.numeric(atol)
  if ( length(rtol) != 1 )
    stop("rtol must (currently) be scalar")
  if ( length(atol) == 1 )
    atol <- rep(atol, n.var)
  else if ( length(atol) != n.var )
    stop("atol must be scalar or of length ", n.var)
  ptr <- .Call("r_make_cvodes", n.var, n.par, derivs, rtol, atol,
               PACKAGE="diversitree")

  function(pars, vars, times) {
    if ( length(pars) != n.par )
      stop("Incorrect parameter length")
    if ( length(vars) != n.var )
      stop("Incorrect variable length")
    if ( length(times) <= 1 )
      stop("Need >= 2 times")
    storage.mode(pars) <- storage.mode(vars) <- storage.mode(times) <-
      "numeric"
    .Call("r_cvodes_set_pars", ptr, pars, PACKAGE="diversitree")
    .Call("r_cvodes_run", ptr, vars, times, PACKAGE="diversitree")
  }
}

cvodes.sens <- function(n.var, n.par, derivs, sens1, rtol, atol,
                        package=NULL) {
  check.cvodes(error=TRUE)

  ## Directly copied from above...
  n.var <- as.integer(n.var)
  n.par <- as.integer(n.par)
  if ( is.null(package) ) package <- ""
  derivs <- getNativeSymbolInfo(derivs, PACKAGE=package)$address
  rtol <- as.numeric(rtol)
  atol <- as.numeric(atol)
  if ( length(rtol) != 1 )
    stop("rtol must (currently) be scalar")
  if ( length(atol) == 1 )
    atol <- rep(atol, n.var)
  else if ( length(atol) != n.var )
    stop("atol must be scalar or of length ", n.var)

  ## New below here:
  sens1 <- getNativeSymbolInfo(sens1, PACKAGE=package)$address
  ptr <- .Call("r_make_cvodes_fwd", n.var, n.par, derivs, sens1,
               rtol, atol, PACKAGE="diversitree")
  n.sens <- n.var * n.par

  function(pars, vars, vars.sens, times, sens=TRUE) {
    ## Directly copied from above:
    if ( length(pars) != n.par )
      stop("Incorrect parameter length")
    if ( length(vars) != n.var )
      stop("Incorrect variable length")
    if ( length(times) <= 1 )
      stop("Need >= 2 times")
    if ( length(vars.sens) != n.sens )
      stop("Incorrect sens length")
    storage.mode(pars) <- storage.mode(vars) <- storage.mode(times) <-
      storage.mode(vars.sens) <- "numeric"
    .Call("r_cvodes_set_pars", ptr, pars, PACKAGE="diversitree")
    ## if (sens) ... else .Call("r_cvodes_run", ptr, vars, times)# ?
    .Call("r_cvodes_fwd_run", ptr, vars, vars.sens, times,
          PACKAGE="diversitree")
  }
}

## See check.fftC()
check.cvodes <- function(error=TRUE) {
  ok <- is.loaded("r_make_cvodes", "diversitree")
  if ( error && !ok )
    stop("diversitree built without CVODES support")
  ok
}
