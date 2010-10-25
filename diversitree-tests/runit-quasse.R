## Tests for QuaSSE likelihood calculations.
test.quasse <- function() {
  ## Build tree:
  lambda <- function(x) sigmoid.x(x, 0.1, 0.2,  0, 2.5)
  mu <- function(x) constant.x(x, 0.03)
  char <- make.brownian.with.drift(0, 0.025)

  set.seed(1)
  phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0,
                     single.lineage=FALSE, verbose=TRUE)

  pars <- c(.1, .2, 0, 2.5, .03, 0, .01)

  sd <- 1/200
  control.C.1 <- list(dt.max=1/200)
  control.C.2 <- c(control.C.1, tips.combined=TRUE)
  control.M.1 <- list(method="mol")
  control.R.1 <- list(dt.max=1/200, method="fftR")

  lik.C.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                         control.C.1)
  lik.C.2 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                         control.C.2)
  lik.M.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                         control.M.1)
  lik.R.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                         control.R.1)

  ll.C.1 <- lik.C.1(pars)
  ll.C.2 <- lik.C.2(pars)
  ll.M.1 <- lik.M.1(pars)
  ll.R.1 <- lik.R.1(pars)

  checkEquals(ll.C.1, -52.90819)
  checkEquals(ll.C.1, ll.C.2)
  checkEquals(ll.C.1, ll.M.1, tolerance=0.0007)
  checkEquals(ll.C.1, ll.R.1)
}
