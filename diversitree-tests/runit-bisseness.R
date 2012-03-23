library(diversitree)
library(RUnit)

test.bisseness <- function() {
  ## First we simulat a 50 species tree, assuming cladogenetic shifts in 
  ## the trait (i.e., the trait only changes at speciation).
  ## Red is state '1', black is state '0', and we let red lineages
  ## speciate at twice the rate of black lineages.
  ## The simulation starts in state 0.
  set.seed(3)
  pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
  phy <- tree.bisseness(pars, max.taxa=50, x0=0)

  set.seed(1)
  pars2 <- pars + runif(10, 0, .1)

  ## This builds the likelihood of the data according to BiSSEness:
  lik1 <- make.bisseness(phy, phy$tip.state)
  lik2 <- make.bisseness(phy, phy$tip.state, control=list(backend="cvodes"))
  ## e.g., the likelihood of the true parameters is:

  checkEquals(lik1(pars), -174.79539045802767)
  checkEquals(lik2(pars), -174.79539045802767, tol=1e-7)

  checkEquals(lik1(pars2), -178.26442717965469)
  checkEquals(lik2(pars2), -178.26442717965469, tol=1e-7)

  set.seed(1)
  pars2 <- pars + runif(10, 0, .1)
  checkEquals(lik1(pars2, root=ROOT.FLAT),
              -178.27220997143513)
  checkEquals(lik1(pars2, root=ROOT.OBS),
              -178.26442717965469)
  checkEquals(lik1(pars2, root=ROOT.EQUI),
              -178.27903183100543)
  checkEquals(lik1(pars2, root=ROOT.FLAT, condition.surv=FALSE),
              -181.72354033281960)
  checkEquals(lik1(pars2, root=ROOT.OBS, condition.surv=FALSE),
              -181.57759200998927)
  checkEquals(lik1(pars2, root=ROOT.EQUI, condition.surv=FALSE),
              -181.83678256438199)

  checkEquals(lik2(pars2, root=ROOT.FLAT),
              -178.27220997143513, tol=1e-7)
  checkEquals(lik2(pars2, root=ROOT.OBS),
              -178.26442717965469, tol=1e-7)
  checkEquals(lik2(pars2, root=ROOT.EQUI),
              -178.27903183100543, tol=1e-7)
  checkEquals(lik2(pars2, root=ROOT.FLAT, condition.surv=FALSE),
              -181.72354033281960, tol=1e-7)
  checkEquals(lik2(pars2, root=ROOT.OBS, condition.surv=FALSE),
              -181.57759200998927, tol=1e-7)
  checkEquals(lik2(pars2, root=ROOT.EQUI, condition.surv=FALSE),
              -181.83678256438199, tol=1e-7)

  ## ## ML search:  First we make hueristic guess at a starting point, based
  ## ## on the constant-rate birth-death model assuming anagenesis (uses
  ## ## \link{make.bd}).
  ## startp <- starting.point.bisse(phy)

  ## ## We then take the total amount of anagenetic change expected across
  ## ## the tree and assign half of this change to anagenesis and half to
  ## ## cladogenetic change at the nodes as a heuristic starting point:
  ## t <- branching.times(phy)
  ## tryq <- 1/2 * startp[["q01"]] * sum(t)/length(t)
  ## p <- c(startp[1:4], startp[5:6]/2, p0c=tryq, p0a=0.5, p1c=tryq, p1a=0.5)

  ## ## Start an ML search from this point.  This takes some time (~12s)
  ## fit <- find.mle(lik, p, method="subplex")
  ## logLik(fit) # -174.0104

  ## ## Compare the fit to a constrained model that only allows the trait
  ## ## to change along a lineage (anagenesis).  This takes some time (~12s)
  ## lik.no.clado <- constrain(lik, p0c ~ 0, p1c ~ 0)
  ## fit.no.clado <- find.mle(lik.no.clado,p[argnames(lik.no.clado)])
  ## logLik(fit.no.clado) # -174.0577

  ## ## This is consistent with what BiSSE finds:
  ## likB <- make.bisse(phy, phy$tip.state)
  ## fitB <- find.mle(likB, startp, method="subplex")
  ## logLik(fitB) # -174.0576

  ## With only this 50-species tree, there is no statistical support
  ## for the more complicated BiSSE-ness model that allows cladogenesis:
  ## anova(fit, no.clado=fit.no.clado)
  ## Note that anova() performs a likelihood ratio test here.

  ## If the above is repeated with max.taxa=250, BiSSE-ness rejects the
  ## constrained model in favor of one that allows cladogenetic change.

  ## Unresolved tip clade: Here we collapse one clade in the 50 species
  ## tree (involving sister species sp70 and sp71) and illustrate the use
  ## of BiSSEness with unresolved tip clades.
  slimphy <- drop.tip(phy,c("sp71"))
  states <- slimphy$tip.state[slimphy$tip.label]
  states["sp70"] <- NA
  unresolved <- data.frame(tip.label=c("sp70"), Nc=2, n0=2, n1=0)

  ## This builds the likelihood of the data according to BiSSEness:
  ## SW: warning suppression here because of the warning about
  ## unresolved clades not being tested extensively.
  lik.unresolved <-
    suppressWarnings(make.bisseness(slimphy, states, unresolved))
  ## e.g., the likelihood of the true parameters is:
  checkEquals(lik.unresolved(pars), -174.65748563924021)
  checkEquals(lik.unresolved(pars2), -178.18221333874183)
}
