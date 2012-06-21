test.bd <- function() {
  library(diversitree)
  library(RUnit)

  set.seed(1)
  phy <- trees(c(.1, .03), "bd", max.taxa=25)[[1]]
  lik <- make.bd(phy)

  lik1 <- make.bd(phy, control=list(method="ode"))
  lik2 <- make.bd(phy, control=list(method="ode", method="cvodes"))
  lik3 <- make.bd(phy, control=list(method="ode", method="CVODES"))

  p <- c(.1, .03)
  checkEquals(lik(p), -22.50266563415211)
  checkEquals(lik1(p), lik(p), tolerance=2e-7)
  checkEquals(lik2(p), lik(p), tolerance=2e-7)
  checkEquals(lik3(p), lik(p), tolerance=2e-7)

  set.seed(1)
  p <- c(.1, .03) + runif(2, 0, .2)
  checkEquals(lik(p),  -23.59180370809460)
  checkEquals(lik1(p), lik(p), tolerance=2e-7)
  checkEquals(lik2(p), lik(p), tolerance=2e-7)
  checkEquals(lik3(p), lik(p), tolerance=2e-7)

  checkEquals(lik(p, FALSE),  -27.26771051981412)
  checkEquals(lik1(p, FALSE), lik(p, FALSE), tolerance=2e-7)
  checkEquals(lik2(p, FALSE), lik(p, FALSE), tolerance=2e-7)
  checkEquals(lik3(p, FALSE), lik(p, FALSE), tolerance=2e-7)

  fit.ape <- birthdeath(phy)
  p.ape <- rev(fit.ape$para) # zero extinction..

  fit.dt <- find.mle(lik, p, method="optim", lower=0)
  checkEquals(coef(fit.dt), p.ape, tolerance=1e-5, check.attr=FALSE)
  checkEquals(-fit.ape$dev/2, fit.dt$lnLik)

  lik.s <- make.bd(phy, sampling.f=.5)
  lik1.s <- make.bd(phy, sampling.f=.5, control=list(method="ode"))
  lik2.s <- make.bd(phy, sampling.f=.5,
                    control=list(method="ode", method="cvodes"))
  lik3.s <- make.bd(phy, sampling.f=.5,
                    control=list(method="ode", method="CVODES"))
  checkEquals(lik.s(p), -23.81317463410290)
  checkEquals(lik1.s(p), lik.s(p), tolerance=2e-7)
  checkEquals(lik2.s(p), lik.s(p), tolerance=2e-7)
  checkEquals(lik3.s(p), lik.s(p), tolerance=2e-7)
  checkEquals(lik1.s(p, FALSE), lik.s(p, FALSE), tolerance=2e-7)
  checkEquals(lik2.s(p, FALSE), lik.s(p, FALSE), tolerance=2e-7)
  checkEquals(lik3.s(p, FALSE), lik.s(p, FALSE), tolerance=2e-7)  
}  

test.bd.t <- function() {
  library(diversitree)
  library(RUnit)

  set.seed(1)
  pars <- c(.1, .03)
  phy <- trees(pars, "bd", max.taxa=25)[[1]]

  ## Next, make three different likelihood functions: a "normal" one that
  ## uses the direct birth-death calculation, an "ode" based one (that
  ## uses numerical integration to compute the likelihood, and is
  ## therefore not exact), and one that is time-varying, but that the
  ## time-dependent functions are constant.t().
  lik.nee <- make.bd(phy)
  lik.ode <- make.bd(phy, control=list(method="ode"))
  lik.t <- make.bd.t.old(phy, list(constant.t, constant.t))

  checkEquals(lik.nee(pars), -22.50266563415211)
  
  ## ODE-based likelihood calculations are correct to about 1e-6.
  checkEquals(lik.ode(pars), lik.nee(pars), tol=2e-7)

  ## The ODE calculation agrees exactly with the time-varying (but
  ## constant) calculation.
  checkIdentical(lik.ode(pars), lik.t(pars))

  ## Next, make a real case, where speciation is a linear function of
  ## time.
  lik.t2 <- make.bd.t.old(phy, list(linear.t, constant.t))

  ## Confirm that this agrees with the previous calculations when the
  ## slope is zero
  pars2 <- c(pars[1], 0, pars[2])
  checkIdentical(lik.t2(pars2),  lik.t(pars))

  set.seed(1)
  pars3 <- pars2 + runif(length(pars2), 0, .2)
  checkEquals(lik.t2(pars3), -105.43147936963197)
}

test.bd.split <- function() {
  set.seed(1)
  pars <- c(.1, .03)
  phy <- trees(pars, "bd", max.taxa=30)[[1]]

  ## Construct the plain likelihood function as a benchmark:
  lik <- make.bd(phy)
  lik(pars) # -21.74554

  ## Split this phylogeny at three points: nd11, nd13 and nd26
  nodes <- c("nd11", "nd13", "nd26")

  ## This is the index in ape's node indexing:
  nodes.i <- match(nodes, phy$node.label) + length(phy$tip.label)

  ## To make a split likelihood function, pass the node locations and times in:
  lik.s <- make.bd.split(phy, nodes)

  ## The parameters must be a list of the same length as the number of
  ## partitions.  Partition '1' is the root partition, and partition i is
  ## the partition rooted at the node[i-1]
  pars4 <- rep(pars, 4)
  names(pars4) <- argnames(lik.s)

  ## Run the likelihod calculation:
  lik.s(pars4) # -21.74554

  ## These are basically identical (to acceptable tolerance)
  checkEquals(lik.s(pars4), lik(pars))

  set.seed(1)
  pars2 <- pars + runif(2, 0, .1)
  checkEquals(lik.s(rep(pars2, 4)), lik(pars2))

  set.seed(1)
  pars4.2 <- pars4 + runif(length(pars4), 0, .1)
  checkEquals(lik.s(pars4.2), -23.48722743120750)

  ## You can use the labelled nodes rather than indices:
  lik.s2 <- make.bd.split(phy, nodes)
  checkIdentical(lik.s(pars4), lik.s2(pars4))

  ## All the usual ML/MCMC functions work as before:
  fit <- find.mle(lik.s, pars4)
  checkEquals(fit$lnLik, -19.65871005552911)
}
