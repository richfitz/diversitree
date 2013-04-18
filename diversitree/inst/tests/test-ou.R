library(diversitree)
library(testthat)
suppressMessages(library(geiger))

no.stdout <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  eval.parent(substitute(expr))
}

context("Ornstein-Uhlenbeck")

data(bird.orders)
set.seed(1)
x <- structure(rnorm(length(bird.orders$tip.label)),
               names=bird.orders$tip.label)

## With the VCV approach
## Note here that the value of theta has no effect on the
## likelihood.  This is a bit surprising, actually.
lik.vcv <- make.ou(bird.orders, x, control=list(method="vcv"))
fit1 <- find.mle(lik.vcv, c(.1, .1, .1))

lik.pruning <- make.ou(bird.orders, x, control=list(method="pruning"))
fit2 <- find.mle(lik.pruning, c(.1, .1, .1))

lik.pruning.C <- make.ou(bird.orders, x,
                         control=list(method="pruning", backend="C"))
fit3 <- find.mle(lik.pruning.C, c(.1, .1, .1))

lik.pruning.C(coef(fit2))

expect_that(fit1, equals(fit2))
## expect_that(fit1, equals(fit3))
expect_that(fit1$lnLik, equals(fit3$lnLik))
expect_that(coef(fit1)[1:2], equals(coef(fit3)[1:2]))

fit4 <- no.stdout(fitContinuous(bird.orders, x, model="OU"))
expect_that(fit4$opt$lnL, equals(fit1$lnLik))
## These are quite different, but the likelihood there suggests a
## ridge...
## expect_that(c(fit4$opt$sigsq,fit4$opt$alpha),
##             equals(fit1$par[1:2], check.attributes=FALSE))
