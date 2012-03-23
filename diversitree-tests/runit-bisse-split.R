test.bisse.split <- function() {
  library(diversitree)
  library(RUnit)

  ## First, simulate the tree:
  pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  set.seed(546)
  phy <- tree.bisse(pars, max.taxa=30, x0=0)

  lik.b <- make.bisse(phy, phy$tip.state)
  lik.b(pars) # -93.62479
  checkEquals(lik.b(pars), -93.62479172896994)

  nodes <- c("nd15", "nd18", "nd26")
  nodes.i <- match(nodes, phy$node.label) + length(phy$tip.label)
  split.t <- Inf # TODO: Drop
  lik.s <- make.bisse.split(phy, phy$tip.state, nodes.i, split.t)

  pars.s <- rep(pars, 4)
     
  ## Run the likelihod calculation:
  checkEquals(lik.s(pars.s), lik.b(pars))

  ## Check the node identity...
  lik.s2 <- make.bisse.split(phy, phy$tip.state, nodes, split.t)
  identical(lik.s(pars.s), lik.s2(pars.s))

  set.seed(1)
  pars2 <- pars + runif(6, 0, .1)
  pars2.s <- rep(pars2, 4)
  checkEquals(lik.b(pars2), -94.16272144707750)
  checkEquals(lik.s(pars2.s), lik.b(pars2))

  pars3.s <- pars + runif(length(pars.s), 0, .1)
  checkEquals(lik.s(pars3.s), -94.18230685862810)
  identical(lik.s(pars3.s), lik.s2(pars3.s))

  unresolved <-
    data.frame(tip.label=c("sp12", "sp32", "sp9", "sp22", "sp11"),
               Nc=c(2,5,3,2,5), n0=c(1, 4, 3, 2, 4), n1=c(1, 1, 0, 0, 1))

  ## Plain BiSSE with unresolved clades:
  lik.u.b <- make.bisse(phy, phy$tip.state, unresolved=unresolved)
  checkEquals(lik.u.b(pars), -139.36879978895746)

  ## Split BiSSE with unresolved clades:
  lik.u.s <- make.bisse.split(phy, phy$tip.state, nodes, split.t,
                              unresolved=unresolved)
  checkEquals(lik.u.s(pars.s),  lik.u.b(pars))
  checkEquals(lik.u.s(pars2.s), lik.u.b(pars2))
  checkEquals(lik.u.s(pars3.s), -127.78056192052314)
}
