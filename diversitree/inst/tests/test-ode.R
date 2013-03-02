## Testing for the ODE bits.
library(diversitree)
library(testthat)

## Here are the BiSSE equations.  This was the first model implemented
## in diversitree, so it makes sense to start here.
derivs <- function(t, y, pars) {
  E0 <- y[1]
  E1 <- y[2]
  D0 <- y[3]
  D1 <- y[4]

  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  c(-(mu0 + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1,
    -(mu1 + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0,
    -(mu0 + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1,
    -(mu1 + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0)
}

derivs.for.desolve <- function(f)
  function(...) list(f(...))
  
## Initial conditions corresponding to a tip in state 0:
y <- c(0, 0, 1, 0)

## Vector of parameters
pars <- c(.1, .2, .03, .04, .01, .02)

## Vector of times to report at
tt <- seq(0, 30, length.out=101)

## Control list for the ODE calculations.
check.control.ode <- diversitree:::check.control.ode
control.deSolve <- check.control.ode(list(backend="deSolve"))
control.gslode  <- check.control.ode(list(backend="gslode"))

## Run the calculations with deSolve's lsoda (no clever tricks)
res.ref <- lsoda(y, tt, derivs.for.desolve(derivs), pars,
                 atol=control.deSolve$tol, rtol=control.deSolve$tol)
## Convert this to the format that we expect (drop time column and
## transpose data).
res.ref <- unname(t(res.ref[,-1]))

## Now build an "info" list:
info <- list(name="bisse", np=6, ny=4, idx.d=3:4,
             argnames=c("la0", "la1", "mu0", "mu1", "q01", "q10"))
info.deSolve.R <- c(info, list(derivs=derivs.for.desolve(derivs)))
info.gslode.R <- c(info, list(derivs=derivs))

## And build the ODE
make.ode <- diversitree:::make.ode
ode.deSolve <- make.ode(info.deSolve.R, control.deSolve)

res.deSolve <- ode.deSolve(y, tt, pars)

expect_that(res.ref[,-1,drop=FALSE], equals(res.deSolve))
expect_that(res.ref[,-1,drop=FALSE], is_identical_to(res.deSolve))

ode.gslode <- make.ode(info.gslode.R, control.gslode)

res.gslode <- ode.gslode(y, tt, pars)
expect_that(res.ref[,-1,drop=FALSE], equals(res.gslode))

## With compiled derivatives:
ode.deSolve.C <- make.ode(info, control.deSolve)
res.deSolve.C <- ode.deSolve.C(y, tt, pars)
expect_that(res.ref[,-1,drop=FALSE], equals(res.deSolve.C))

ode.gslode.C <- make.ode(info, control.gslode)
res.gslode.C <- ode.gslode.C(y, tt, pars)
expect_that(res.ref[,-1,drop=FALSE], equals(res.gslode.C))

## Now, try a time dependent model.
## Really, need to sort out doing the time machine augmenting
## somewhere else...

functions <- rep(c("linear.t", "constant.t"), c(2, 4))
names(functions) <- info$argnames

## Drop t.range from this?  Or do something with it?
info$tm <- diversitree:::make.time.machine2(functions, c(0, 10))
info$time.varying <- TRUE
info$argnames <- attr(info$tm, "argnames")
info$np <- length(info$argnames)

ode.gslode.t <- make.ode(info, control.gslode)

pars.t <- c(pars[1], 0, pars[2], 0, pars[3:6])

res.gslode.t <- ode.gslode.t(y, tt, pars.t)
expect_that(res.ref[,-1,drop=FALSE], equals(res.gslode.t))

## Cool - that does appear to work.  Now let's add in some time
## dependence properly:
pars.t2 <- pars.t
pars.t2[c(2, 4)] <- 0.01
res.gslode.t2 <- ode.gslode.t(y, tt, pars.t2)

## We'll need a BiSSE equation to compare this with:
## Here are the BiSSE equations.  This was the first model implemented
## in diversitree, so it makes sense to start here.
derivs.t <- function(t, y, pars) {
  E0 <- y[1]
  E1 <- y[2]
  D0 <- y[3]
  D1 <- y[4]

  lambda0 <- pars[1] + pars[2] * t
  lambda1 <- pars[3] + pars[4] * t
  mu0 <- pars[5]
  mu1 <- pars[6]
  q01 <- pars[7]
  q10 <- pars[8]

  c(-(mu0 + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1,
    -(mu1 + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0,
    -(mu0 + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1,
    -(mu1 + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0)
}

res.ref.t <- lsoda(y, tt, derivs.for.desolve(derivs.t), pars.t2,
                   atol=control.deSolve$tol, rtol=control.deSolve$tol)
## Convert this to the format that we expect (drop time column and
## transpose data).
res.ref.t <- unname(t(res.ref.t[,-1]))

expect_that(res.ref.t[,-1,drop=FALSE], equals(res.gslode.t2))


