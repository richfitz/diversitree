protect <- function(f) {
  function(..., fail.value=NULL, finite=TRUE) {
    if ( is.null(fail.value) )
      f(...)
    else {
      ret <- try(f(...), silent = TRUE)
      failed <- (inherits(ret, "try-error") ||
                 (finite && (is.na(ret) || !is.finite(ret))))
      if ( failed )
        fail.value
      else
        ret
    }
  }
}

invert <- function(f) function(...) -f(...)

big.brother <- function(f, interval=1) {
  .x.eval <- list()
  .y.eval <- list()
  function(x, ...) {
    i <- length(.x.eval) + 1
    if ( i %% interval == 0 )
      cat(sprintf("[%s]", paste(formatC(x, 5), collapse=", ")))
    else
      cat(".")
    .x.eval[[i]] <<- x
    .y.eval[[i]] <<- ans <- f(x, ...)
    if ( i %% interval == 0 )
      cat(sprintf("\t -> %2.5f\n", ans))
    ans
  }
}

## This is updated for deSolve version 1.5, but now includes a "safe"
## argument that tries to avoid doing any silly business.
make.ode <- function(func, dllname, initfunc, ny, safe=FALSE) {
  if (!is.character(func)) 
    stop("`func' must be a character vector")
  if (!is.character(dllname))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")

  if ( safe ) {
    function(y, times, parms, rtol, atol)
      lsoda(y, times, parms, rtol=rtol, atol=atol,
            initfunc=initfunc, func=func, dllname=dllname)
  } else {
    ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE=dllname)$address
    Func <- getNativeSymbolInfo(func, PACKAGE=dllname)$address
    JacFunc <- NULL

    maxordn <- 12
    maxords <- 5
    lrw <- as.integer(max(20 + ny * (maxordn + 1) + 3 * ny,
                          20 + ny * (maxords + 1) + 3 * ny + ny^2 + 2))
    liw <- as.integer(20 + ny)

    iwork <- vector("integer", 20)
    iwork[1] <- as.integer(1)    # banddown
    iwork[2] <- as.integer(1)    # bandup
    iwork[6] <- as.integer(5000) # maxsteps

    rwork <- vector("double", 20) # TODO: 5, 6 do nothing.
    rwork[5] <- 0                     # hini
    rwork[6] <- 10                    # hmax (consider 0)
    rwork[7] <- 0                     # hmin

    INTZERO <- as.integer(0)
    INTONE <- as.integer(1)
    INTTWO <- as.integer(2)

    ## "Forcings" are new to deSolve version 1.5; this is the required
    ## argument where they are not used (I am not using them)
    flist <- list(fmat=0, tmat=0, imat=0, ModelForc=NULL)  
    
    function(y, times, parms, rtol, atol) {
      if (!is.numeric(y)) 
        stop("`y' must be numeric")
      if (!is.numeric(times)) 
        stop("`times' must be numeric")
      storage.mode(y) <- storage.mode(times) <- "double"
      .Call("call_lsoda", y, times, Func, parms, rtol, atol,
            NULL,    # rho: environment
            NULL,    # tcrit: critical times
            JacFunc,
            ModelInit,
            INTZERO, # verbose (false)
            INTONE,  # itask (no idea)
            rwork,
            iwork,
            INTTWO,  # jt: Jacobian type (fullint)
            INTZERO, # Nglobal (no global variables)
            lrw,     # size of workspace (real)
            liw,     # size of workspace (int)
            INTONE,  # 'IN' (no idea)
            NULL,    # ?
            INTZERO, # ?
            0,       # rpar: no extra real parameters
            INTZERO, # ipar: no extra integer parameters
            INTZERO, # ?
            flist,   # New for deSolve 1.5
            PACKAGE="deSolve")
    }
  }
}

quadratic.roots <- function(a, b, c)
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)

## TODO: Do the tips immediately after the loop; this should be really
## easy, but I think that I need to include the edge matrix to make it
## work.  Should be something like this:
##   ans[tips] <- ans[match(tips, edge[,2])]
ancestors <- function(parent, order) {
  ans <- vector("list", max(order))
  for ( i in rev(order[-length(order)]) )
    ans[[i]] <- c(parent[i], ans[[parent[i]]])
  ans
}
