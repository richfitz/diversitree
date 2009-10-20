## This is not specific to BiSSE.
make.ode <- function(func, dllname, initfunc, ny) {
  if (!is.character(func)) 
    stop("`func' must be a character vector")
  if (!is.character(dllname))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")

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
          PACKAGE="deSolve")
  }
}

## Here is the R version:
bisse.eqs <- function(t, x, pars) {
  E0 <- x[1]
  E1 <- x[2]
  D0 <- x[3]
  D1 <- x[4]

  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  list(c(-(mu0 + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1,
         -(mu1 + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0,
         -(mu0 + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1,
         -(mu1 + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0))
}

## The Jacobian is not used, but is here for reference:
make.jacobian <- function(t, x, pars) {
  E0 <- x[1]
  E1 <- x[2]
  D0 <- x[3]
  D1 <- x[4]

  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  jac <- matrix(0, 4, 4)
  jac[1,1] <- -(mu0 + q01 + lambda0) + 2 * lambda0 * E0
  jac[1,2] <- q01

  jac[2,1] <- q10
  jac[2,2] <- -(mu1 + q10 + lambda1) + 2 * lambda1 * E1

  jac[3,1] <- 2 * D0 * lambda0
  jac[3,3] <- -(mu0 + q01 + lambda0) + 2 * lambda0 * E0
  jac[3,4] <- q01

  jac[4,2] <- 2 * D1 * lambda1
  jac[4,3] <- q01
  jac[4,4] <- -(mu1 + q10 + lambda1) + 2 * lambda1 * E1
  jac
}

solve.R <- function(y, t, pars)
  t(ode(y, c(0, t), bisse.eqs, pars))[-1,]

## Default tolerances
RTOL <- 1e-8
ATOL <- 1e-8
solve.C <- function(y, t, pars, rtol=RTOL, atol=ATOL)
  bisse.ode(y, c(0, t), pars, rtol=rtol, atol=atol)[,-1]

solve <- solve.C

