## Known to work on 1.3-1.6
make.ode.deSolve <- function(info, control) {
  model  <- info$name.ode
  np     <- info$np  
  ny     <- info$ny
  dll    <- info$dll  
  safe   <- control$safe
  unsafe <- control$unsafe
  atol   <- rtol <- as.numeric(control$tol)

  initfunc <- sprintf("initmod_%s", model)
  derivs <- sprintf("derivs_%s", model)

  if ( safe ) {
    function(vars, times, pars) {
      ret <- t(lsoda(vars, times, pars, rtol=rtol, atol=atol,
                     initfunc=initfunc, func=derivs,
                     dllname=dll)[-1,-1,drop=FALSE])
      dimnames(ret) <- NULL
      ret
    }
  } else {
    ## Temporary fix so that I can work on the cluster.  This will be
    ## removed and DESCRIPTION updated to require R 2.12.0 or greater.
    if ( getRversion() >= package_version("2.12.0") )
      vers <-  packageVersion("deSolve")
    else
      vers <- package_version(packageDescription("deSolve",
                                                 fields="Version"))
    max.deSolve <- package_version("1.10-3")
    if ( !unsafe && vers > max.deSolve ) {
      str <- paste("diversitree is not known to work with deSolve > ",
                   max.deSolve, "\n\tfalling back on slow version")
      warning(str)
      control$safe <- TRUE
      return(make.ode.deSolve(info, control))
    }
    
    initfunc  <- getNativeSymbolInfo(initfunc, PACKAGE=dll)$address
    derivs    <- getNativeSymbolInfo(derivs,   PACKAGE=dll)$address
    jacfunc <- NULL

    maxordn <- 12
    maxords <- 5
    lrw <- as.integer(max(20 + ny * (maxordn + 1) + 3 * ny,
                          20 + ny * (maxords + 1) + 3 * ny + ny^2 + 2))
    liw <- as.integer(20 + ny)

    iwork <- vector("integer", 20)
    iwork[1] <- as.integer(1)    # banddown
    iwork[2] <- as.integer(1)    # bandup
    iwork[6] <- as.integer(5000) # maxsteps

    rwork <- vector("double", 20)
    rwork[5] <- 0                     # hini
    rwork[6] <- 10                    # hmax (consider 0)
    rwork[7] <- 0                     # hmin

    INTZERO <- 0L
    INTONE <- 1L
    INTTWO <- 2L

    flist <- list(fmat=0, tmat=0, imat=0, ModelForc=NULL)
    elag <- list(islag=0L)

    sol <- function(vars, times, pars) {
      if ( length(vars) != ny )
        stop("Incorrect variable length")
      if ( length(times) <= 1 )
        stop("Need >= 2 times")
      storage.mode(vars) <- storage.mode(times) <- "numeric"

      ret <- 
        .Call("call_lsoda", vars, times, derivs, pars,
              rtol, atol,
              NULL,      # rho: environment
              NULL,      # tcrit: critical times
              jacfunc, 
              initfunc,
              NULL,      # eventfunc [New in 1.6]
              INTZERO,   # verbose (false)
              INTONE,    # itask
              rwork,
              iwork,
              INTTWO,    # jT: Jacobian type (fullint)
              INTZERO,   # nOut (no global variables)
              lrw,       # size of workspace (real)
              liw,       # size of workspace (int)
              INTONE,    # Solver
              NULL,      # rootfunc
              INTZERO,   # nRoot
              0,         # rpar: no extra real parameters
              INTZERO,   # ipar: no extra integer parameters
              INTZERO,   # Type
              flist,     # [New in 1.5]
              list(),    # elist [New in 1.6]
              elag,      # [New in 1.7]
              PACKAGE="deSolve")
      if ( max(abs(ret[1,] - times)) > 1e-6 )
        stop("Integration error: integration stopped prematurely")
      ret[-1,-1,drop=FALSE]
    }
  }
}

## For testing new models.  Selected by make.ode() when info$derivs is
## a function
make.ode.R <- function(info, control) {
  derivs <- info$derivs
  rtol <- atol <- control$tol
  if ( !is.function(derivs) )
    stop("info$derivs must be a function")

  function(vars, times, pars) {
    ret <- t(lsoda(vars, times, derivs, pars,
                   rtol=rtol, atol=atol)[-1,-1,drop=FALSE])
    dimnames(ret) <- NULL
    ret
  }
}
