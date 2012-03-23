test.asr.bisse <- function() {
  library(diversitree)
  library(RUnit)
  
  pars <- c(.1, .2, .03, .06, .01, .02)
  set.seed(3)
  phy <- trees(pars, "bisse", max.taxa=50, max.t=Inf, x0=0)[[1]]

  lik <- make.bisse(phy, phy$tip.state)
  p <- c(0.084, 0.211, 0.061, 0, 0.004, 0.024)
  st <- asr.marginal(lik, p)

  checkEquals(var(st[1,]), 0.14547319392772)
  checkEquals(st[1,3], 0.58689689938945)

  set.seed(1)
  p2 <- p + runif(length(p), 0, .1)
  st2 <- asr.marginal(lik, p2)
  checkEquals(var(st2[1,]), 0.08703485887271)
  checkEquals(st2[1,3], 0.06446214107037)

  lik.m <- make.musse(phy, phy$tip.state+1, 2)
  checkEquals(asr.marginal(lik.m, p2), st2)

  lik.m <- make.mk2(phy, phy$tip.state)
  fit.m <- find.mle(lik.m, pars[5:6], method="subplex")
  st.m <- asr.marginal(lik.m, coef(fit.m))

  st.id <- asr.marginal(lik, c(.1, .1, .03, .03, coef(fit.m)))
  checkEquals(st.id, st.m, tolerance=1e-7)  

  set.seed(1)
  st.id2 <- asr.marginal(lik, c(rep(runif(2, 0, .1), each=2), coef(fit.m)))
  checkEquals(st.id2, st.m, tolerance=1e-7)  

  lik2 <- make.bisse(phy, phy$tip.state, control=list(backend="cvodes"))
  lik3 <- make.bisse(phy, phy$tip.state, control=list(backend="CVODES"))

  st2 <- asr.marginal(lik2, p)
  checkEquals(st, st2, tol=5e-7)

  st3 <- asr.marginal(lik3, p)
  checkEquals(st, st3, tol=5e-7)
}

test.asr.musse <- function() {
  library(diversitree)
  library(RUnit)

  pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
            .03, .045, .06, # mu 1, 2, 3
            .05, 0,         # q12, q13
            .05, .05,       # q21, q23
            0,   .05)       # q31, q32
  set.seed(2)
  phy <- tree.musse(pars, 30, x0=1)

  states <- phy$tip.state  
  lik <- make.musse(phy, states, 3)

  set.seed(1)
  p <- pars + runif(length(pars), 0, .2)
  
  checkEquals(lik(p), -118.95287844068999)
  st <- asr.marginal(lik, p)

  checkEquals(var(st[1,]), 0.05216825252108)
  checkEquals(st[1,3], 0.05180212059737)

  lik.m <- make.mkn(phy, states, 3)
  p2 <- p
  p2[1:3] <- .1
  p2[4:6] <- .03
  p2.m <- p2[-(1:6)]
  
  checkEquals(lik(p2), -119.15072625137572)
  checkEquals(lik.m(p2.m), -28.64702621988359)

  st1 <- asr.marginal(lik, p2)
  st2 <- asr.marginal(lik.m, p2.m)
  checkEquals(st1, st2, tol=5e-7)
}
