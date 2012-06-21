## Check that the ancestral state estimation under time-dependent
## models makes sense.

#### 1. BiSSE (easier)
test.asr.t.bisse <- function() {
  library(diversitree)
  library(RUnit)

  ## Start with a simple tree evolved under a BiSSE with all rates
  ## asymmetric:
  pars <- c(.1, .2, .03, .06, .01, .02)
  set.seed(3)
  phy <- trees(pars, "bisse", max.taxa=50, max.t=Inf, x0=0)[[1]]

  ## BiSSE ancestral state reconstructions under the ML model
  lik.b <- make.bisse(phy, phy$tip.state)
  lik.m <- make.musse(phy, phy$tip.state+1, 2)
  ## fit <- find.mle(lik, pars, method="subplex")
  ## p.b.ml <- dput(coef(fit))

  p.b <- pars
  p.b.ml <- structure(c(0.0840286576638916, 0.211042245168045,
                        0.0607003387416787, 3.34756696282173e-06,
                        0.00412941736207307, 0.0236232224532126),
                      .Names=c("lambda0", "lambda1", "mu0", "mu1", "q01",
                        "q10"))

  st.b <- asr.marginal(lik.b, p.b)
  st.m <- asr.marginal(lik.m, p.b)
  checkEquals(st.b, st.m)

  lik.b.t <- make.bisse.t(phy, phy$tip.state,
                          rep(c("linear.t", "constant.t"), c(2, 4)))
  lik.b.t2 <- make.bisse.t(phy, phy$tip.state,
                           rep(c("linear.t", "constant.t"), c(2, 4)),
                           control=list(backend="cvodes"))
  lik.b.t3 <- make.bisse.t(phy, phy$tip.state,
                           rep(c("linear.t", "constant.t"), c(2, 4)),
                           control=list(backend="CVODES"))
  lik.m.t <- make.musse.t(phy, phy$tip.state+1, 2,
                          rep(c("linear.t", "constant.t"), c(2, 4)))
  lik.m.t2 <- make.musse.t(phy, phy$tip.state+1, 2,
                           rep(c("linear.t", "constant.t"), c(2, 4)),
                           control=list(backend="cvodes"))
  lik.m.t3 <- make.musse.t(phy, phy$tip.state+1, 2,
                           rep(c("linear.t", "constant.t"), c(2, 4)),
                           control=list(backend="CVODES"))

  ## Time-independent parameters for time models
  p.b.t <- c(rbind(p.b[1:2], 0), p.b[-(1:2)])
  names(p.b.t) <- argnames(lik.b.t)

  ## Time-dependent parameters for time models
  p.b.tt <- p.b.t
  p.b.tt[sprintf("lambda%d.m", 0:1)] <- .002

  ## Sanity check -- make sure that the likelihoods agree.
  ll.t <- lik.b.t(p.b.t)
  checkEquals(ll.t, lik.b.t2(p.b.t))
  checkEquals(ll.t, lik.b.t3(p.b.t), tolerance=1e-7)
  checkEquals(ll.t, lik.m.t(p.b.t))
  checkEquals(ll.t, lik.m.t2(p.b.t))
  checkEquals(ll.t, lik.m.t3(p.b.t), tolerance=1e-7)

  ll.tt <- lik.b.t(p.b.tt)
  checkEquals(ll.tt, lik.b.t2(p.b.tt), tolerance=4e-7)
  checkEquals(ll.tt, lik.b.t3(p.b.tt), tolerance=4e-7)
  checkEquals(ll.tt, lik.m.t(p.b.tt))
  checkEquals(ll.tt, lik.m.t2(p.b.tt), tolerance=4e-7)
  checkEquals(ll.tt, lik.m.t3(p.b.tt), tolerance=4e-7)

  ## Check the ancestral states:
  st.b.t <- asr.marginal(lik.b.t, p.b.t)
  st.b.t2 <- asr.marginal(lik.b.t2, p.b.t)
  st.b.t3 <- asr.marginal(lik.b.t3, p.b.t)
  st.m.t <- asr.marginal(lik.m.t, p.b.t)
  st.m.t2 <- asr.marginal(lik.m.t2, p.b.t)
  st.m.t3 <- asr.marginal(lik.m.t3, p.b.t)

  checkIdentical(st.b, st.b.t)
  checkIdentical(st.m, st.m.t)

  checkEquals(st.b, st.b.t2, tol=5e-7)
  checkEquals(st.b, st.b.t3, tol=5e-7)
  checkEquals(st.m, st.m.t2, tol=5e-7)
  checkEquals(st.m, st.m.t3, tol=5e-7)

  st.b.tt <- asr.marginal(lik.b.t, p.b.tt)
  st.b.tt2 <- asr.marginal(lik.b.t2, p.b.tt)
  st.b.tt3 <- asr.marginal(lik.b.t3, p.b.tt)
  st.m.tt <- asr.marginal(lik.m.t, p.b.tt)
  st.m.tt2 <- asr.marginal(lik.m.t2, p.b.tt)
  st.m.tt3 <- asr.marginal(lik.m.t3, p.b.tt)

  checkEquals(st.b.tt, st.b.tt2, tolerance=1e-6)
  checkEquals(st.b.tt, st.b.tt3, tolerance=1e-6)
  checkEquals(st.b.tt, st.m.tt, tolerance=1e-6)
  checkEquals(st.b.tt, st.m.tt2, tolerance=1e-6)
  checkEquals(st.b.tt, st.m.tt3, tolerance=1e-6)
}

test.asr.t.musse <- function() {
#### 2. MuSSE (harder)
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
  table(states)

  lik.m <- make.musse(phy, states, 3)
  lik.t0 <- make.musse.t.old(phy, states, 3,
                             rep(c(linear.t, constant.t), c(3, 9)))
  lik.t1 <- make.musse.t(phy, states, 3,
                         rep(c("linear.t", "constant.t"), c(3, 9)))
  lik.t2 <- make.musse.t(phy, states, 3,
                         rep(c("linear.t", "constant.t"), c(3, 9)),
                         control=list(backend="cvodes"))
  lik.t3 <- make.musse.t(phy, states, 3,
                         rep(c("linear.t", "constant.t"), c(3, 9)),
                         control=list(backend="CVODES"))

  p <- starting.point.musse(phy, 3)
  p.t <- c(rbind(p[1:3], 0), p[-(1:3)])
  names(p.t) <- argnames(lik.t1)

  ## Likelihoods agree.
  checkEquals(lik.m(p), lik.t0(p.t))
  checkEquals(lik.m(p), lik.t1(p.t))
  checkEquals(lik.m(p), lik.t2(p.t), tol=1e-7)
  checkEquals(lik.m(p), lik.t3(p.t), tol=1e-7)

  ## Ancestral states:
  st.m <- asr.marginal(lik.m, p)
  ## st.t0 <- asr.marginal(lik.t0, p.t) # no longer works.
  st.t1 <- asr.marginal(lik.t1, p.t)
  st.t2 <- asr.marginal(lik.t2, p.t)
  st.t3 <- asr.marginal(lik.t3, p.t)

  checkIdentical(st.m, st.t1)
  checkEquals(st.m, st.t2, tol=1e-6)
  checkEquals(st.m, st.t3, tol=1e-6)

  p.tt <- p.t
  p.tt[sprintf("lambda%d.m", 1:3)] <- c(.02, .015, .01)

  st.tt1 <- asr.marginal(lik.t1, p.tt)
  st.tt2 <- asr.marginal(lik.t2, p.tt)
  st.tt3 <- asr.marginal(lik.t3, p.tt)

  checkEquals(rowSums(st.tt1),
              c(11.4584593797639, 13.3218167093813, 4.21972391085483))
  checkEquals(rowSums(st.tt2), rowSums(st.tt1), tol=2e-6)
  checkEquals(rowSums(st.tt3), rowSums(st.tt1), tol=2e-6)
}
