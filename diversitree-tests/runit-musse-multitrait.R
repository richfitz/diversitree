test.musse.multitrait <- function() {
  library(diversitree)
  library(RUnit)
  
  tr <- musse.multitrait.translate(2)
  pars <- c(.10, .15, .20, .25, # lambda 00, 10, 01, 11
            .03, .03, .03, .03, # mu 00, 10, 01, 11
            .05, .05, .0,       # q00.10, q00.01, q00.11
            .05, .0,  .05,      # q10.00, q10.01, q10.11
            .05, .0,  .05,      # q01.00, q01.10, q01.11
            .0,  .05, .05)      # q11.00, q11.10, q11.01
  set.seed(2)
  phy <- tree.musse(pars, 60, x0=1)

  states <- expand.grid(A=0:1, B=0:1)[phy$tip.state,]
  rownames(states) <- phy$tip.label

  ## Note that transition from the original MuSSE basis to this basis is
  ## only possible in general when depth=n.trait and allow.multistep=TRUE
  ## (as only this generates a square matrix that is invertible).
  ## However, when it is possible to express the set of parameters in the
  ## new basis (as it is above), this can be done through a pseudoinverse
  ## (here, a left inverse).
  pars2 <- drop(solve(t(tr) %*% tr) %*% t(tr) %*% pars)

  ## Going from our new basis to the original MuSSE parameters is always
  ## straightforward.  This is done automatically in the likelihood
  ## function.
  checkEquals(c(tr %*% pars2), pars, check.attributes=FALSE)

  ## This shows that the two traits act additively on speciation rate
  ## (lambdaAB is zero), that there is no effect of any trait on
  ## extinction (the only nonzero mu parameter is mu0) and transition
  ## rates for one trait are unaffected by other traits (the only nonzero
  ## q parameters are the qXij.0 parameters; qXij.Y parameters are all
  ## zero).

  ## Here is our new MuSSE function parametrised as a multi-trait
  ## function:
  lik1 <- make.musse.multitrait(phy, states)
  lik2 <- make.musse.multitrait(phy, states,
                                control=list(backend="cvodes"))
  lik3 <- make.musse.multitrait(phy, states,
                                control=list(backend="CVODES"))

  ## Basic MuSSE function for comparison
  lika <- make.musse(phy, phy$tip.state, 4)
  likb <- make.musse(phy, phy$tip.state, 4,
                     control=list(backend="cvodes"))
  likc <- make.musse(phy, phy$tip.state, 4,
                     control=list(backend="CVODES"))

  checkEquals(lik1(pars2), lika(pars))
  checkEquals(lik1(pars2, root=ROOT.FLAT),
              lika(pars,  root=ROOT.FLAT))
  checkEquals(lik1(pars2, root=ROOT.OBS),
              lika(pars,  root=ROOT.OBS))
  checkEquals(lik1(pars2, root=ROOT.FLAT, condition.surv=TRUE),
              lika(pars,  root=ROOT.FLAT, condition.surv=TRUE))
  checkEquals(lik1(pars2, root=ROOT.OBS,  condition.surv=TRUE),
              lika(pars,  root=ROOT.OBS,  condition.surv=TRUE))

  checkEquals(lik2(pars2), likb(pars))
  checkEquals(lik2(pars2, root=ROOT.FLAT),
              likb(pars,  root=ROOT.FLAT))
  checkEquals(lik2(pars2, root=ROOT.OBS),
              likb(pars,  root=ROOT.OBS))
  checkEquals(lik2(pars2, root=ROOT.FLAT, condition.surv=TRUE),
              likb(pars,  root=ROOT.FLAT, condition.surv=TRUE))
  checkEquals(lik2(pars2, root=ROOT.OBS,  condition.surv=TRUE),
              likb(pars,  root=ROOT.OBS,  condition.surv=TRUE))

  checkEquals(lik3(pars2), likc(pars))
  checkEquals(lik3(pars2, root=ROOT.FLAT),
              likc(pars,  root=ROOT.FLAT))
  checkEquals(lik3(pars2, root=ROOT.OBS),
              likc(pars,  root=ROOT.OBS))
  checkEquals(lik3(pars2, root=ROOT.FLAT, condition.surv=TRUE),
              likc(pars,  root=ROOT.FLAT, condition.surv=TRUE))
  checkEquals(lik3(pars2, root=ROOT.OBS,  condition.surv=TRUE),
              likc(pars,  root=ROOT.OBS,  condition.surv=TRUE))
  
  pars.t <- pars2 + runif(length(pars2), 0, .1)
  pars.m <- c(tr %*% pars.t)

  checkEquals(lik1(pars.t), lika(pars.m))
  checkEquals(lik1(pars.t, root=ROOT.FLAT),
              lika(pars.m, root=ROOT.FLAT))
  checkEquals(lik1(pars.t, root=ROOT.OBS),
              lika(pars.m, root=ROOT.OBS))
  checkEquals(lik1(pars.t, root=ROOT.FLAT, condition.surv=TRUE),
              lika(pars.m, root=ROOT.FLAT, condition.surv=TRUE))
  checkEquals(lik1(pars.t, root=ROOT.OBS,  condition.surv=TRUE),
              lika(pars.m, root=ROOT.OBS,  condition.surv=TRUE))

  checkEquals(lik2(pars.t), likb(pars.m))
  checkEquals(lik2(pars.t, root=ROOT.FLAT),
              likb(pars.m, root=ROOT.FLAT))
  checkEquals(lik2(pars.t, root=ROOT.OBS),
              likb(pars.m, root=ROOT.OBS))
  checkEquals(lik2(pars.t, root=ROOT.FLAT, condition.surv=TRUE),
              likb(pars.m, root=ROOT.FLAT, condition.surv=TRUE))
  checkEquals(lik2(pars.t, root=ROOT.OBS,  condition.surv=TRUE),
              likb(pars.m, root=ROOT.OBS,  condition.surv=TRUE))

  checkEquals(lik3(pars.t), likc(pars.m))
  checkEquals(lik3(pars.t, root=ROOT.FLAT),
              likc(pars.m, root=ROOT.FLAT))
  checkEquals(lik3(pars.t, root=ROOT.OBS),
              likc(pars.m, root=ROOT.OBS))
  checkEquals(lik3(pars.t, root=ROOT.FLAT, condition.surv=TRUE),
              likc(pars.m, root=ROOT.FLAT, condition.surv=TRUE))
  checkEquals(lik3(pars.t, root=ROOT.OBS,  condition.surv=TRUE),
              likc(pars.m, root=ROOT.OBS,  condition.surv=TRUE))
}
