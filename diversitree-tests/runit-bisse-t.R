test.bisse.t <- function() {
  library(diversitree)
  library(RUnit)
  set.seed(4)
  pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  phy <- tree.bisse(pars, max.t=30, x0=0)

  lik.b <- make.bisse(phy, phy$tip.state)  
  lik.t <- make.bisse.t.old(phy, phy$tip.state,
                            rep(list(sigmoid.t, constant.t), c(2, 4)))

  ## First, evaluate the functions with no time effect and check that they
  ## are the same as the base BiSSE model
  t <- max(branching.times(phy))/2
  p.t <- c(pars[1], pars[1], t, .5,
           pars[2], pars[2], t, .5,
           pars[3:6])

  checkEquals(lik.b(pars), -159.709955958895989)
  checkEquals(lik.b(pars), lik.t(p.t))

  set.seed(1)
  pars2 <- pars + runif(6, 0, .2)
  p2.t <- c(pars2[1], pars2[1], t, .5,
            pars2[2], pars2[2], t, .5,
            pars2[3:6])
  
  checkEquals(lik.b(pars2), -189.60809396398534)
  checkEquals(lik.b(pars2), lik.t(p2.t))

  p3.t <- p2.t
  p3.t[2] <- p3.t[2] + .2
  checkEquals(lik.t(p3.t), -190.58110120320532)

  ## Now, using the new interface:
  lik.n <- make.bisse.t(phy, phy$tip.state,
                        rep(c("sigmoid.t", "constant.t"), c(2, 4)))
  checkEquals(lik.n(p.t), lik.t(p.t))
  checkEquals(lik.n(p2.t), lik.t(p2.t))
  checkEquals(lik.n(p3.t), lik.t(p3.t))

  lik.n2 <- make.bisse.t(phy, phy$tip.state,
                         rep(c("sigmoid.t", "constant.t"), c(2, 4)),
                         control=list(backend="cvodes"))
  checkEquals(lik.n2(p.t), lik.t(p.t), tolerance=2e-7)
  checkEquals(lik.n2(p2.t), lik.t(p2.t), tolerance=2e-7)
  checkEquals(lik.n2(p3.t), lik.t(p3.t), tolerance=2e-7)
  
  lik.N <- make.bisse.t(phy, phy$tip.state,
                        rep(c("sigmoid.t", "constant.t"), c(2, 4)),
                        control=list(backend="CVODES"))

  checkEquals(lik.n2(p.t), lik.N(p.t), tolerance=2e-7)
  checkEquals(lik.n2(p2.t), lik.N(p2.t), tolerance=2e-7)
  checkEquals(lik.n2(p3.t), lik.N(p3.t), tolerance=2e-7)
}

test.bisse.td <- function() {
  library(diversitree)
  library(RUnit)
  set.seed(4)
  pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  phy <- tree.bisse(pars, max.t=30, x0=0)

  lik.b <- make.bisse(phy, phy$tip.state)
  lik.td <- make.bisse.td(phy, phy$tip.state, 2)
  
  ## First, evaluate the functions with no time effect and check that they
  ## are the same as the base BiSSE model
  t <- max(branching.times(phy))/2
  p.td <- c(t, pars, pars)

  checkEquals(lik.b(pars), -159.709955958895989)
  checkEquals(lik.b(pars), lik.td(p.td))

  set.seed(1)
  pars2 <- pars + runif(6, 0, .2)
  p2.td <- c(t, pars2, pars2)  
  
  checkEquals(lik.b(pars2), -189.60809396398534)
  checkEquals(lik.b(pars2), lik.td(p2.td))

  set.seed(1)
  p3.td <- p2.td + runif(length(p2.td), 0, .2)
  checkEquals(lik.td(p3.td), -183.12543259220922)

  lik.t <- make.bisse.t.old(phy, phy$tip.state, rep(list(stepf.t), 6))
  t <- p3.td[1]
  p <- matrix(p3.td[-1], ncol=2)
  p3.t <- c(p[1,], t,
            p[2,], t,
            p[3,], t,
            p[4,], t,
            p[5,], t,
            p[6,], t)
  names(p3.t) <- argnames(lik.t)
  names(p3.td) <- argnames(lik.td)
  checkEquals(lik.t(p3.t), -183.12543259220922)
}

test.musse.t <- function() {
  library(diversitree)
  library(RUnit)
  pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
            .03, .045, .06, # mu 1, 2, 3
            .05, 0,         # q12, q13
            .05, .05,       # q21, q23
            0,   .05)       # q31, q32
  set.seed(2)
  phy <- tree.musse(pars, 30, x0=1)

  ## For comparison, make a plain MuSSE likelihood function
  lik.m <- make.musse(phy, phy$tip.state, 3)

  lik.t <- make.musse.t.old(phy, phy$tip.state, 3,
                            rep(list(sigmoid.t, constant.t), c(3, 9)))

  ## First, evaluate the functions with no time effect and check that they
  ## are the same as the base MuSSE model
  t <- max(branching.times(phy))/2
  p.t <- c(pars[1], pars[1], t, .5,
           pars[2], pars[2], t, .5,
           pars[3], pars[3], t, .5,
           pars[4:12])

  checkEquals(lik.m(pars), -110.836439594910701)
  checkEquals(lik.m(pars), lik.t(p.t))

  set.seed(1)
  pars2 <- pars + runif(length(pars), 0, .2)
  p2.t <- c(pars2[1], pars2[1], t, .5,
            pars2[2], pars2[2], t, .5,
            pars2[3], pars2[3], t, .5,
            pars2[4:12])
  
  checkEquals(lik.m(pars2), -118.952878440689986)
  checkEquals(lik.m(pars2), lik.t(p2.t))

  p3.t <- p2.t
  p3.t[2] <- p3.t[2] + .2
  checkEquals(lik.t(p3.t), -118.752042795013892)
}

test.musse.td <- function() {
  library(diversitree)
  library(RUnit)
  pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
  set.seed(2)
  phy <- tree.musse(pars, 30, x0=1)

  ## For comparison, make a plain MuSSE likelihood function
  lik.m <- make.musse(phy, phy$tip.state, 3)

  ## Create the time-dependent likelihood function.  The final argument
  ## here is the number of 'epochs' that are allowed.  Two epochs is one
  ## switch point.
  lik.td <- make.musse.td(phy, phy$tip.state, 3, 2)

  argnames(lik.td)

  t <- max(branching.times(phy))/2  
  p.td <- c(t, pars, pars)
  names(p.td) <- argnames(lik.td)

  checkEquals(lik.m(pars), -110.836439594910701)
  checkEquals(lik.m(pars), lik.td(p.td))

  set.seed(1)
  pars2 <- pars + runif(length(pars), 0, .2)
  p2.td <- c(t, pars2, pars2)
  checkEquals(lik.m(pars2), -118.952878440689986)
  checkEquals(lik.m(pars2), lik.td(p2.td))

  set.seed(1)
  p3.td <- p2.td + runif(length(p2.td), 0, .2)
  checkEquals(lik.td(p3.td), -122.52408021988634)
}
