test.musse.t <- function() {
  library(diversitree)
  library(RUnit)
  ## 1: BiSSE equivalence
  pars <- c(.1, .2, .03, .04, 0.05, 0.1)
  set.seed(2)
  phy <- tree.musse(pars, 20, x0=1)

  ## Show that the likelihood functions give the same answers.  Ignore the
  ## warning when creating the MuSSE function.
  lik.m <- make.musse(phy, phy$tip.state, 2)

  lik.b <- make.bisse.t(phy, phy$tip.state-1,
                        rep(c("linear.t", "constant.t"), c(2, 4)))
  lik.t <- make.musse.t(phy, phy$tip.state, 2,
                        rep(c("linear.t", "constant.t"), c(2, 4)))

  pars.t <- c(pars[1], 0, pars[2], 0, pars[-(1:2)])

  ## Everything checks out when time-independent.
  checkEquals(lik.b(pars.t), lik.m(pars), tolerance=1e-7)
  checkEquals(lik.m(pars), lik.t(pars.t))

  pars.t2 <- pars.t
  set.seed(1)
  pars.t2[c(2,4)] <- runif(2)*0.01
  checkEquals(lik.b(pars.t2), lik.t(pars.t2))

  ## You can get the expected argument order for any number of states
  ## this way (sorry - clunky).  The help file also lists the order.
  diversitree:::default.argnames.musse(3)

  ## Here are the parameters:
  pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
            .03, .045, .06, # mu 1, 2, 3
            .05, 0,         # q12, q13
            .05, .05,       # q21, q23
            0,   .05)       # q31, q32

  set.seed(2)
  phy <- tree.musse(pars, 30, x0=1)

  ## The states are numbered 1:3, rather than 0:1 in bisse.
  states <- phy$tip.state
  table(states)

  lik.m <- make.musse(phy, states, 3)
  lik.t <- make.musse.t(phy, states, 3,
                        rep(c("linear.t", "constant.t"), c(3, 9)))
  lik.T <- make.musse.t(phy, states, 3,
                        rep(c("linear.t", "constant.t"), c(3, 9)),
                        control=list(backend="CVODES"))
  lik.old <- make.musse.t.old(phy, states, 3,
                              rep(c(linear.t, constant.t), c(3, 9)))

  p <- starting.point.musse(phy, 3)
  p.t <- c(rbind(p[1:3], 0), p[-(1:3)])
  names(p.t) <- argnames(lik.t)

  checkEquals(lik.m(p), lik.t(p.t))

  p.t2 <- p.t
  set.seed(1)
  p.t2[c(2,4,6)] <- runif(3)*0.01
  checkEquals(lik.t(p.t2), -111.9320762)
  checkEquals(lik.T(p.t2), -111.9320762, tolerance=1e-7)
  checkEquals(lik.old(p.t2), -111.9320762)

  ## Next, all in with a variable model
  lik.t2 <- make.musse.t(phy, states, 3, "linear.t")
  p.t3 <- c(rbind(pars, 0))
  names(p.t3) <- argnames(lik.t2)
  checkEquals(lik.t2(p.t3), lik.m(pars))

  set.seed(1)
  p.t4 <- c(rbind(pars, runif(length(pars))*0.01))
  names(p.t4) <- argnames(lik.t2)
  checkEquals(lik.t2(p.t4), -113.780614532)
}
