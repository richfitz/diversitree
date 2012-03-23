test.musse <- function() {
  library(diversitree)
  library(RUnit)
  ## 1: BiSSE equivalence
  pars <- c(.1, .2, .03, .04, 0.05, 0.1)
  set.seed(2)
  phy <- tree.musse(pars, 20, x0=1)

  ## Show that the likelihood functions give the same answers.  Ignore the
  ## warning when creating the MuSSE function.
  lik.b <- make.bisse(phy, phy$tip.state-1)
  lik.m <- make.musse(phy, phy$tip.state, 2)
  checkEquals(lik.b(pars), lik.m(pars), tolerance=1e-7)

  ## Notice that default argument names are different between BiSSE and
  ## MuSSE, but that the order is the same.
  argnames(lik.b) # BiSSE: 0/1
  argnames(lik.m) # MuSSE: 1/2

  ## 2: A 3-state example where movement is only allowed between
  ## neighbouring states (1 <-> 2 <-> 3), and where speciation and
  ## extinction rates increase moving from 1 -> 2 -> 3:

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

  ## 2: Likelihood
  ## Making a likelihood function is basically identical to bisse.  The
  ## third argument needs to be the number of states.  In a future
  ## version this will probably be max(states), but there are some
  ## pitfalls about this that I am still worried about.
  lik <- make.musse(phy, states, 3)

  ## Here are the arguments.  Even with three states, this is getting
  ## ridiculous.
  argnames(lik)

  ## Start with a fully constrained model, but still enforcing stepwise
  ## changes (disallowing 1 <-> 3 shifts)
  lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1,
                        mu2 ~ mu1, mu3 ~ mu1,
                        q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)

  ## TODO: Error here of all places...
  p <- starting.point.musse(phy, 3)
  ## Should be OK once this is done.
  fit.base <- find.mle(lik.base, p[argnames(lik.base)])
  checkEquals(fit.base$lnLik, -110.634754)

  ## Now, allow the speciation rates to vary:
  lik.lambda <- constrain(lik, mu2 ~ mu1, mu3 ~ mu1,
                          q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
  fit.lambda <- find.mle(lik.lambda, p[argnames(lik.lambda)])
  checkEquals(fit.lambda$lnLik, -110.2117391)
}

test.musse.split <- function() {
  ## First, simulate the tree:
  set.seed(2)
  pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
            .03, .045, .06, # mu 1, 2, 3
            .05, 0,         # q12, q13
            .05, .05,       # q21, q23
            0,   .05)       # q31, q32
  phy <- tree.musse(pars, 30, x0=1)

  ## Here is a plain MuSSE function for later comparison:
  lik.m <- make.musse(phy, phy$tip.state, 3)
  lik.m(pars) # -110.8364
  checkEquals(lik.m(pars), -110.8364396)

  ## Split this phylogeny at three points: nd16 and nd25, splitting it
  ## into three chunks
  nodes <- c("nd16", "nd25")

  ## To make a split BiSSE function, pass the node locations and times
  ## in.  Here, we'll use 'Inf' as the split time to mimick MEDUSA's
  ## behaviour of placing the split at the base of the branch subtended by
  ## a node.
  lik.s <- make.musse.split(phy, phy$tip.state, 3, nodes, split.t=Inf)

  ## The parameters must be a list of the same length as the number of
  ## partitions.  Partition '1' is the root partition, and partition i is
  ## the partition rooted at the node[i-1]:
  argnames(lik.s)

  ## Because we have two nodes, there are three sets of parameters.
  ## Replicate the original list to get a starting point for the analysis:
  pars.s <- rep(pars, 3)
  names(pars.s) <- argnames(lik.s)

  checkEquals(lik.s(pars.s), -110.8364396)

  ## This is basically identical (to acceptable tolerance) to the plain
  ## MuSSE version:
  checkEquals(lik.s(pars.s), lik.m(pars))

  set.seed(1)
  pars.s2 <- pars.s + runif(36, 0, .1)
  checkEquals(lik.s(pars.s2), -114.292614038858)

  ## The resulting likelihood function can be used in ML analyses with
  ## find.mle.  However, because of the large number of parameters, this
  ## may take some time (especially with as few species as there are in
  ## this tree - getting convergence in a reasonable number of iterations
  ## is difficult).

  ## Bayesian analysis also works, using the mcmc function.  Given the
  ## large number of parameters, priors will be essential, as there will
  ## be no signal for several parameters.  Here, I am using an exponential
  ## distribution with a mean of twice the state-independent
  ## diversification rate.
}
