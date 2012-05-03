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

test.ou <- function() {
  library(diversitree)
  library(RUnit)
  data(bird.orders)
  set.seed(1)
  x <- structure(rnorm(length(bird.orders$tip.label)),
                 names=bird.orders$tip.label)

  ## With the VCV approach
  lik.vcv <- make.ou(bird.orders, x, control=list(method="vcv"))
  system.time(fit1 <- find.mle(lik.vcv, c(.1, .1, .1)))

  lik.pruning <- make.ou(bird.orders, x, control=list(method="pruning"))
  system.time(fit2 <- find.mle(lik.pruning, c(.1, .1, .1)))

  lik.pruning.C <- make.ou(bird.orders, x,
                           control=list(method="pruning", backend="C"))
  system.time(fit3 <- find.mle(lik.pruning.C, c(.1, .1, .1)))

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

