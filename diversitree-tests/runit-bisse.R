## Low level BiSSE tests
test.bisse <- function() {
  pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  set.seed(4)
  phy <- tree.bisse(pars, max.t=30, x0=0)

  ## Here is the plain likelihood calculation.
  lik <- make.bisse(phy, phy$tip.state)
  lik(pars) # -159.71

  ## Here, generate an unresolved description that is the same as the
  ## individual tips (i.e., single species, in observed state)
  unresolved <- data.frame(tip.label=phy$tip.label, Nc=1,
                           n0=1-phy$tip.state, n1=phy$tip.state)
  states <- phy$tip.state
  states[] <- NA

  ## This gives the same warning that you are seeing:
  lik2 <- make.bisse(phy, states, unresolved=unresolved)

  ## But the likelihood calculations are through both functions are
  ## identical to a reasonable tolerance (each individual numerical
  ## integration is only accurate to a relative tolerance of about 1e-8)
  lik2(pars) - lik(pars)

  checkEquals(lik2(pars), lik(pars), 1e-6)
}
