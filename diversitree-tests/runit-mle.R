## Really excercise the MLE code
test.mle <- function() {
  ## Start with a simple 2 parameter model from the BD model
  pars <- c(0.1, 0.03)
  set.seed(2)
  phy <- tree.bd(pars, max.taxa=60)
  lik <- make.bd(phy)

  ## First, exercise some basic options:
  ans.nlm.1 <- find.mle(lik, pars, method="nlm")
  ans.nlm.2 <- find.mle(lik, pars, method="nlm", verbose=40)
  ans.nlm.3 <- find.mle(lik, pars, method="nlm", fail.value=-500)
  ans.nlm.4 <- find.mle(lik, pars, method="nlm",
                        control=list(fail.penalty=500))

  checkIdentical(ans.nlm.1, ans.nlm.2)
  checkEquals(ans.nlm.1[1:2], ans.nlm.3[1:2], tolerance=1e-2)
  checkIdentical(!identical(ans.nlm.1, ans.nlm.4), TRUE)
  checkEquals(ans.nlm.1[1:2], ans.nlm.4[1:2], tolerance=0.07)

  ## Next, excercise the different methods:
  ans.optim <- suppressWarnings(find.mle(lik, pars, method="optim"))
  ans.subplex <- find.mle(lik, pars, method="subplex")
  ans.nlminb <- find.mle(lik, pars, method="nlminb")
  ans.nlm <- find.mle(lik, pars, method="nlm")
  ans.minqa <- find.mle(lik, pars, method="minqa")

  ## Different optim methods:
  ans.optim.nm <- find.mle(lik, pars, method="optim",
                           control=list(optim.method="Nelder-Mead"))
  ans.optim.bfgs <- find.mle(lik, pars, method="optim",
                             control=list(optim.method="BFGS"))
  ans.optim.cg <- find.mle(lik, pars, method="optim",
                           control=list(optim.method="CG"))
  ans.optim.lbfgsb <-
    suppressWarnings(find.mle(lik, pars, method="optim",
                              control=list(optim.method="L-BFGS-B")))
  ans.optim.sann <- find.mle(lik, pars, method="optim",
                             control=list(optim.method="SANN"))

  ## Different minqa methods:
  ans.minqa.n <- find.mle(lik, pars, method="minqa",
                          control=list(minqa.method="newuoa"))
  ans.minqa.b <- find.mle(lik, pars, method="minqa",
                          control=list(minqa.method="bobyqa"))
  ans.minqa.u <- find.mle(lik, pars, method="minqa",
                          control=list(minqa.method="uobyqa"))

  ## "Odd" methods, involving integers.

  ## Here is a "likelihood" function based on the Rosenbrock banana
  ## function.  We will evaluate this function only at integer values of
  ## x.
  lik <- function(x, x2, as.integer=TRUE) {
    x1 <- if ( as.integer ) diversitree:::check.integer(x) else x
    -(100*(x2-x1*x1)^2+(1-x1)^2)
  }

  fit <- find.mle(lik, -11, x2=-33, method="int1d")
  checkIdentical(fit$par, 0)

  ## And a few options here, too:
  fit.1 <- find.mle(lik, -11, x2=-33, method="int1d", upper=11)
  fit.2 <- find.mle(lik, -11, x2=-33, method="int1d", lower=-11)
  fit.3 <- find.mle(lik, -11, x2=-33, method="int1d",
                    control=list(interval=c(-11, 1)))
  checkIdentical(fit.1$par, 0)
  checkIdentical(fit.2$par, 0)
  checkIdentical(fit.3$par, 0)

  ## Now, the mixed method: the 4d rosenbrock function with some integer
  ## axes.
  rosen.multi <- function(x) {
    xx <- matrix(x, 2)
    x1 <- xx[1,]
    x2 <- xx[2,]
    -sum(100*(x2-x1*x1)^2+(1-x1)^2)
  }

  ## Let's have the third argument be an integer.
  rosen.mixed <- function(x) {
    diversitree:::check.integer(x[3])
    rosen.multi(x)
  }

  ## This works well.
  p <- rep(-11, 4)
  fit.c <- find.mle(rosen.multi, p, method="subplex")

  ## However, this will fail to fit, as the integer axis is highly
  ## correlated with one real axis (in fact, this confuses subplex, and
  ## the other rosenbrock function also does not complete).
  fit.m <- find.mle(rosen.mixed, p, method="mixed",
                    control=list(is.integer=3))
  fit.m2 <- find.mle(rosen.mixed, fit.m$par, method="mixed",
                     control=list(is.integer=3))
  checkIdentical(fit.m2$par[3], -11)
}
