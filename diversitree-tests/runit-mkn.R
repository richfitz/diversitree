test.mkn <- function() {
  library(diversitree)
  library(RUnit)
  ## Simulate a tree and character distribution.  This is on a birth-death
  ## tree, with high rates of character evolution and an asymmetry in the
  ## character transition rates.
  pars <- c(.1, .1, .03, .03, .1, .2)
  set.seed(3)
  phy <- trees(pars, "bisse", max.taxa=25, max.t=Inf, x0=0)[[1]]

  ## Maximum likelihood parameter estimation:
  p <- c(.1, .1) # initial parameter guess
  lik <- make.mk2(phy, phy$tip.state)
  fit.mk2 <- find.mle(lik, p)
  checkEquals(fit.mk2$lnLik, -10.90569519216388)

  lik.mkn <- make.mkn(phy, phy$tip.state + 1, 2)
  fit.mkn <- find.mle(lik.mkn, p)
  checkEquals(fit.mkn$lnLik, -10.90569519216388)

  model <- matrix(c(0, 2, 1, 0), 2)
  fit.ape <- ace(phy$tip.state, phy, "discrete", model=model, ip=p)

  fit.mk2 <- find.mle(lik, p, root=ROOT.GIVEN, root.p=c(1,1))

  checkEquals(fit.ape[c("rates", "loglik")], fit.mk2[1:2],
              check.attributes=FALSE, tolerance=1e-4)

  lik.ode <- make.mkn(phy, phy$tip.state + 1, 2,
                      control=list(method="ode"))
  fit.ode <- find.mle(lik.ode, p)

  checkEquals(fit.ode[-7], fit.mkn[-7], tolerance=1e-7)
}
