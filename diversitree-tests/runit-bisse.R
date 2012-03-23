## Low level BiSSE tests
test.bisse <- function() {
  library(diversitree)
  library(RUnit)
  pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  set.seed(4)
  phy <- tree.bisse(pars, max.t=30, x0=0)

  ## Here is the plain likelihood calculation.
  lik1 <- make.bisse(phy, phy$tip.state)
  lik2 <- make.bisse(phy, phy$tip.state, control=list(backend="cvodes"))
  lik3 <- make.bisse(phy, phy$tip.state, control=list(backend="CVODES"))

  checkEquals(lik1(pars), -159.709955958896)
  checkEquals(lik2(pars), -159.709955958896, tolerance=1e-7)
  checkEquals(lik3(pars), -159.709955958896, tolerance=1e-7)

  set.seed(1)
  pars2 <- pars + runif(6, 0, .1)
  checkEquals(lik1(pars2, root=ROOT.FLAT),
              -173.403159254299)
  checkEquals(lik1(pars2, root=ROOT.OBS),
              -173.083837291932)
  checkEquals(lik1(pars2, root=ROOT.EQUI),
              -173.574847849880)
  checkEquals(lik1(pars2, root=ROOT.FLAT, condition.surv=FALSE),
              -176.564186725815)
  checkEquals(lik1(pars2, root=ROOT.OBS, condition.surv=FALSE),
              -175.939144822762)
  checkEquals(lik1(pars2, root=ROOT.EQUI, condition.surv=FALSE),
              -176.840650405926)

  checkEquals(lik2(pars2, root=ROOT.FLAT),
              -173.403159254299, tolerance=1e-7)
  checkEquals(lik2(pars2, root=ROOT.OBS),
              -173.083837291932, tolerance=1e-7)
  checkEquals(lik2(pars2, root=ROOT.EQUI),
              -173.574847849880, tolerance=1e-7)
  checkEquals(lik2(pars2, root=ROOT.FLAT, condition.surv=FALSE),
              -176.564186725815, tolerance=1e-7)
  checkEquals(lik2(pars2, root=ROOT.OBS, condition.surv=FALSE),
              -175.939144822762, tolerance=1e-7)
  checkEquals(lik2(pars2, root=ROOT.EQUI, condition.surv=FALSE),
              -176.840650405926, tolerance=1e-7)

  checkEquals(lik3(pars2, root=ROOT.FLAT),
              -173.403159254299, tolerance=1e-7)
  checkEquals(lik3(pars2, root=ROOT.OBS),
              -173.083837291932, tolerance=1e-7)
  checkEquals(lik3(pars2, root=ROOT.EQUI),
              -173.574847849880, tolerance=1e-7)
  checkEquals(lik3(pars2, root=ROOT.FLAT, condition.surv=FALSE),
              -176.564186725815, tolerance=1e-7)
  checkEquals(lik3(pars2, root=ROOT.OBS, condition.surv=FALSE),
              -175.939144822762, tolerance=1e-7)
  checkEquals(lik3(pars2, root=ROOT.EQUI, condition.surv=FALSE),
              -176.840650405926, tolerance=1e-7)
  
  ## Here, generate an unresolved description that is the same as the
  ## individual tips (i.e., single species, in observed state)
  unresolved <- data.frame(tip.label=phy$tip.label, Nc=1,
                           n0=1-phy$tip.state, n1=phy$tip.state)
  states <- phy$tip.state
  states[] <- NA

  lik.u <- make.bisse(phy, states, unresolved=unresolved)

  ## But the likelihood calculations are through both functions are
  ## identical to a reasonable tolerance (each individual numerical
  ## integration is only accurate to a relative tolerance of about 1e-8)
  checkEquals(lik.u(pars), lik1(pars))
  checkEquals(lik.u(pars2), lik1(pars2), tol=3e-6)
}
