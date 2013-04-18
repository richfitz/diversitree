library(diversitree)
library(testthat)
suppressMessages(library(geiger))

no.stdout <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  eval.parent(substitute(expr))
}

context("Brownian motion")

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

## All the same:
expect_that(fit1, equals(fit2))
expect_that(fit1, equals(fit3))

library(geiger)
fit4 <- no.stdout(fitContinuous(bird.orders, x))

expect_that(fit1$lnLik, equals(fit4$opt$lnL))
expect_that(fit1$par, equals(fit4$opt$sigsq, tolerance=1e-4,
                             check.attributes=FALSE))

## Now, with standard errors:
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
fit.g.1 <- no.stdout(fitContinuous(phy, x))
fit.g.2 <- no.stdout(fitContinuous(phy, x, SE=se))

## All agree -- nice.
expect_that(unname(coef(fit.v.1)), equals(3.3765625000))
expect_that(coef(fit.p.1), equals(coef(fit.v.1)))
expect_that(coef(fit.P.1), equals(coef(fit.v.1)))
expect_that(lik.v.1(fit.g.1$opt$sigsq), equals(fit.g.1$opt$lnL))

expect_that(unname(coef(fit.v.2)), equals(0.1237304688))
expect_that(coef(fit.p.2), equals(coef(fit.v.2)))
expect_that(coef(fit.P.2), equals(coef(fit.v.2)))
expect_that(lik.v.2(fit.g.2$opt$sigsq), equals(fit.g.2$opt$lnL))

