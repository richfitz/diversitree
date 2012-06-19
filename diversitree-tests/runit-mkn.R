test.mkn <- function() {
  library(diversitree)
  library(RUnit)
  ## Simulate a tree and character distribution.  This is on a birth-death
  ## tree, with high rates of character evolution and an asymmetry in the
  ## character transition rates.
  pars <- c(.1, .1, .03, .03, .1, .2)
  set.seed(3)
  phy <- trees(pars, "bisse", max.taxa=25, max.t=Inf, x0=0)[[1]]

  ## Maximum likelihood parameter estimation:
  p <- c(.1, .1) # initial parameter guess
  lik <- make.mk2(phy, phy$tip.state)
  fit.mk2 <- find.mle(lik, p)
  checkEquals(fit.mk2$lnLik, -10.90569519216388)

  lik.mkn <- make.mkn(phy, phy$tip.state + 1, 2)
  fit.mkn <- find.mle(lik.mkn, p)
  checkEquals(fit.mkn$lnLik, -10.90569519216388)

  model <- matrix(c(0, 2, 1, 0), 2)
  fit.ape <- ace(phy$tip.state, phy, "discrete", model=model, ip=p)

  fit.mk2 <- find.mle(lik, p, root=ROOT.GIVEN, root.p=c(1,1))

  checkEquals(fit.ape[c("rates", "loglik")], fit.mk2[1:2],
              check.attributes=FALSE, tolerance=1e-4)

  lik.ode <- make.mkn(phy, phy$tip.state + 1, 2,
                      control=list(method="ode"))
  fit.ode <- find.mle(lik.ode, p)

  checkEquals(fit.ode[-7], fit.mkn[-7], tolerance=1e-7)

  ## Yay!  This works perfectly for 50, but not always.  I get error
  ## -42, which apparently I must recover from.  FFS.  Not in the main
  ## expokit code?  Is this one of mine?
  ##
  ## This is not much faster for k=50, but 10x faster for k=100.  Only
  ## 20% of time in R code, so possibly not much speed up from a
  ## CVODES style interface.
  ##
  ## meristic is about the same speed, so a Jacobian might be enough.
  set.seed(3)
  k <- 200
  states <- sim.character(phy, c(k, 3, 3), model="meristic", x0=4)
  lik.ode <- make.mkn(phy, states, k, strict=FALSE,
                      control=list(method="ode"))
  lik.meristic <- make.mkn.meristic(phy, states, k)
  lik.expokit <- make.mkn(phy, states, k, strict=FALSE,
                          control=list(method="ode", backend="expokit"))
  pars <- diversitree:::mkn.meristic.Q(c(3, 3), k)
  pars2 <- pars[row(pars) != col(pars)]
  ll.o <- lik.ode(pars2)
  ll.e <- lik.expokit(pars2)
  checkEquals(ll.o, ll.e, tol=1e-7)
}

test.mkn.meristic <- function() {
  library(diversitree)
  library(RUnit)
  ## Simulate a tree and character distribution.  This is on a birth-death
  ## tree, with high rates of character evolution and an asymmetry in the
  ## character transition rates.
  set.seed(3)
  phy <- tree.yule(.1, max.taxa=100)
  k <- 31
  r <- 3
  states <- sim.character(phy, c(k, r, r), x0=(k+1)/2,
                          model="meristic")

  ## Constrain this:
  p <- rep(0, k * (k-1))
  tr <- expand.grid(from=seq_len(k), to=seq_len(k))
  tr <- tr[tr$from != tr$to,]
  base <- ceiling(log10(k + 0.5))
  fmt <- sprintf("q%%0%dd%%0%dd", base, base)
  par <- sprintf(fmt, tr$from, tr$to)
  i.u <- tr$from - tr$to ==  1
  i.d <- tr$from - tr$to == -1
  i.o <- abs(tr$from - tr$to) > 1
  u <- d <- r
  p[i.d] <- d
  p[i.u] <- u
  
  lik.mkn1 <- make.mkn(phy, states, k, strict=FALSE)
  system.time(ll.mkn1 <- lik.mkn1(p))

  lik.mkn2 <- make.mkn(phy, states, k, strict=FALSE,
                       control=list(method="ode"))
  system.time(ll.mkn2 <- lik.mkn2(p))

  lik.mkn3 <- make.mkn(phy, states, k, strict=FALSE,
                       control=list(method="ode", backend="cvodes"))
  system.time(ll.mkn3 <- lik.mkn3(p))
  
  ## I have idea why this is *slower*, and not much faster than .3?
  lik.mkn4 <- make.mkn(phy, states, k, strict=FALSE,
                       control=list(method="ode", backend="CVODES"))
  system.time(ll.mkn4 <- lik.mkn4(p))
  
  lik.mer1 <- make.mkn.meristic(phy, states, k)
  system.time(ll.mer1 <- lik.mer1(c(r, r)))

  lik.mer2 <- make.mkn.meristic(phy, states, k,
                                control=list(backend="cvodes"))
  system.time(ll.mer2 <- lik.mer2(c(d, u)))

  ## This one is also not much faster.  Looks like LSODA just
  ## slaughters CVODE here.  Probably need to tune the problem.
  lik.mer3 <- make.mkn.meristic(phy, states, k,
                                control=list(backend="CVODES"))
  system.time(ll.mer3 <- lik.mer3(c(r, r)))

  checkEquals(ll.mkn1, ll.mkn2, tol=5e-7)
  checkEquals(ll.mkn1, ll.mkn3, tol=5e-7)
  checkEquals(ll.mkn1, ll.mkn4, tol=5e-7)

  checkEquals(ll.mkn1, ll.mer1, tol=5e-7)
  checkEquals(ll.mkn1, ll.mer2, tol=5e-7)
  checkEquals(ll.mkn1, ll.mer3, tol=5e-7)

  ## ## ASYMMETRIC RATES
  u <- 2 * r
  p[i.d] <- d
  p[i.u] <- u
  
  lik.mkn1 <- make.mkn(phy, states, k, strict=FALSE)
  system.time(ll.mkn1 <- lik.mkn1(p))

  lik.mkn2 <- make.mkn(phy, states, k, strict=FALSE,
                       control=list(method="ode"))
  system.time(ll.mkn2 <- lik.mkn2(p))

  lik.mkn3 <- make.mkn(phy, states, k, strict=FALSE,
                       control=list(method="ode", backend="cvodes"))
  system.time(ll.mkn3 <- lik.mkn3(p))
  
  lik.mkn4 <- make.mkn(phy, states, k, strict=FALSE,
                       control=list(method="ode", backend="CVODES"))
  system.time(ll.mkn4 <- lik.mkn4(p))
  
  lik.mer1 <- make.mkn.meristic(phy, states, k)
  system.time(ll.mer1 <- lik.mer1(c(d, u)))

  lik.mer2 <- make.mkn.meristic(phy, states, k,
                                control=list(backend="cvodes"))
  system.time(ll.mer2 <- lik.mer2(c(d, u)))

  ## This one is also not much faster.  Looks like LSODA just
  ## slaughters CVODE here.  Probably need to tune the problem.
  lik.mer3 <- make.mkn.meristic(phy, states, k,
                                control=list(backend="CVODES"))
  system.time(ll.mer3 <- lik.mer3(c(d, u)))

  checkEquals(ll.mkn1, ll.mkn2, tol=5e-5)
  checkEquals(ll.mkn2, ll.mkn3, tol=5e-5)
  checkEquals(ll.mkn2, ll.mkn4, tol=5e-5)

  checkEquals(ll.mkn2, ll.mer1, tol=5e-5)
  checkEquals(ll.mkn2, ll.mer2, tol=5e-5)
  checkEquals(ll.mkn2, ll.mer3, tol=5e-5)

  q1 <- diversitree:::mkn.Q(p)
  q2 <- diversitree:::mkn.meristic.Q(c(d, u), k)
  checkIdentical(q1, q2)  
}

