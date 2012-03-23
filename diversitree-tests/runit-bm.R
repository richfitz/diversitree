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

  ## With the direct calculations
  lik.direct <- make.bm(bird.orders, x, control=list(method="direct"))
  fit2 <- find.mle(lik.direct, .1)

  ## All the same (need to drop the function from this though)
  checkEquals(fit1[-7], fit2[-7])

  library(geiger)
  fit3 <- fitContinuous(bird.orders, x)
  checkEquals(fit3$Trait1$lnl, fit1$lnLik)
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
  lik.ou <- make.ou(bird.orders, x)
  fit1 <- find.mle(lik.ou, c(.1, .1, .1))

  library(geiger)
  fit3 <- fitContinuous(bird.orders, x, model="OU")
  checkEquals(fit3$Trait1$lnl, fit1$lnLik)
  ## These are quite different, but the likelihood there suggests a
  ## ridge...
  ## checkEquals(c(fit3$Trait1$beta,fit3$Trait1$alpha),
  ##             fit1$par[1:2], check.attributes=FALSE)
}

