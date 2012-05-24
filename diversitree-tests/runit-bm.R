test.bm <- function() {
  library(diversitree)
  library(RUnit)
  data(bird.orders)
  set.seed(1)
  x <- structure(rnorm(length(bird.orders$tip.label)),
                 names=bird.orders$tip.label)

  ## With the VCV approach
  lik.vcv <- make.bm(bird.orders, x)
  lik.vcv(.1)
  fit1 <- find.mle(make.bm(bird.orders, x), .1,
                   control=list(method="vcv"))

  ## With the pruning calculations
  lik.pruning <- make.bm(bird.orders, x, control=list(method="pruning"))
  fit2 <- find.mle(lik.pruning, .1)

  ## With the C-based pruning calculations
  lik.pruning.C <- make.bm(bird.orders, x,
                           control=list(method="pruning", backend="C"))
  fit3 <- find.mle(lik.pruning.C, .1)

  all.equal(fit2, fit3)
  
  ## All the same:
  checkEquals(fit1, fit2)
  checkEquals(fit1, fit3)

  library(geiger)
  fit4 <- fitContinuous(bird.orders, x)
  checkEquals(fit4$Trait1$lnl, fit1$lnLik)
  ## These actually differ very slightly:
  ## checkEquals(fit3$Trait1$beta, fit1$par, check.attributes=FALSE)
}

test.bm.stderr <- function() {
  library(diversitree)
  library(geiger)
  library(RUnit)

  set.seed(1)
  phy <- trees(c(.1, .03), "bd", max.taxa=50)[[1]]
  states <- sim.character(phy, .1)
  states["sp62"] <- -10.5
  se <- .8523

  control.v <- list(method="vcv")
  control.p <- list(method="pruning", backend="R")
  control.P <- list(method="pruning", backend="C")

  lik.v.1 <- make.bm(phy, states, 0,  control=control.v)
  lik.v.2 <- make.bm(phy, states, se, control=control.v)
  lik.p.1 <- make.bm(phy, states, 0,  control=control.p)
  lik.p.2 <- make.bm(phy, states, se, control=control.p)
  lik.P.1 <- make.bm(phy, states, 0,  control=control.P)
  lik.P.2 <- make.bm(phy, states, se, control=control.P)

  fit.v.1 <- find.mle(lik.v.1, .1)
  fit.p.1 <- find.mle(lik.p.1, .1)
  fit.P.1 <- find.mle(lik.P.1, .1)

  fit.v.2 <- find.mle(lik.v.2, .1)
  fit.p.2 <- find.mle(lik.p.2, .1)
  fit.P.2 <- find.mle(lik.P.2, .1)

  x <- states
  attr(x, "node.state") <- NULL
  fit.g.1 <- fitContinuous(phy, x)
  fit.g.2 <- fitContinuous(phy, x, meserr=se)

  ## All agree -- nice.
  checkEquals(unname(coef(fit.v.1)), 3.3765625000)
  checkEquals(coef(fit.p.1), coef(fit.v.1))
  checkEquals(coef(fit.P.1), coef(fit.v.1))
  checkEquals(lik.v.1(fit.g.1$Trait1$beta), fit.g.1$Trait1$lnl)

  checkEquals(unname(coef(fit.v.2)), 0.1237304688)
  checkEquals(coef(fit.p.2), coef(fit.v.2))
  checkEquals(coef(fit.P.2), coef(fit.v.2))
  checkEquals(lik.v.2(fit.g.2$Trait1$beta), fit.g.2$Trait1$lnl)
}

test.ou <- function() {
  library(diversitree)
  library(RUnit)
  data(bird.orders)
  set.seed(1)
  x <- structure(rnorm(length(bird.orders$tip.label)),
                 names=bird.orders$tip.label)

  ## With the VCV approach
  ## Note here that the value of theta has no effect on the
  ## likelihood.  This is a bit surprising, actually.
  lik.vcv <- make.ou(bird.orders, x, control=list(method="vcv"))
  system.time(fit1 <- find.mle(lik.vcv, c(.1, .1, .1)))

  lik.pruning <- make.ou(bird.orders, x, control=list(method="pruning"))
  system.time(fit2 <- find.mle(lik.pruning, c(.1, .1, .1)))

  lik.pruning.C <- make.ou(bird.orders, x,
                           control=list(method="pruning", backend="C"))
  system.time(fit3 <- find.mle(lik.pruning.C, c(.1, .1, .1)))

  lik.pruning.C(coef(fit2))

  all.equal(fit1, fit2)
  all.equal(fit1, fit3)
  
  library(geiger)
  system.time(fit4 <- fitContinuous(bird.orders, x, model="OU"))
  checkEquals(fit4$Trait1$lnl, fit1$lnLik)
  ## These are quite different, but the likelihood there suggests a
  ## ridge...
  ## checkEquals(c(fit3$Trait1$beta,fit3$Trait1$alpha),
  ##             fit1$par[1:2], check.attributes=FALSE)
}
