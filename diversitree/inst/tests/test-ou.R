library(diversitree)
library(testthat)
suppressMessages(library(geiger))

no.stdout <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  eval.parent(substitute(expr))
}

context("Ornstein-Uhlenbeck")

## Simulated tree and traits:
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
states <- sim.character(phy, 1)
se <- 0.1

lik.vcv <- make.ou(phy, states, control=list(method="vcv"))
lik.pru.R <- make.ou(phy, states,
                     control=list(method="pruning", backend="R"))
lik.pru.C <- make.ou(phy, states,
                     control=list(method="pruning", backend="C"))
lik.con <- make.ou(phy, states, control=list(method="contrasts"))

lik.vcv.se <- make.ou(phy, states, se, control=list(method="vcv"))
lik.pru.R.se <- make.ou(phy, states, se,
                     control=list(method="pruning", backend="R"))
lik.pru.C.se <- make.ou(phy, states, se,
                     control=list(method="pruning", backend="C"))
## Not yet supported
expect_that(make.ou(phy, states, se,
                    control=list(method="contrasts")),
            throws_error())

## First, a simple test on a value that is known.  Obviously if the
## tree simulators change this will break!  But all four will then
## also break.
test_that("Likelihood calculations agree on known case", {
  ## Start with BM; set alpha to zero:
  pars <- c(1, 0, 0)
  ll <- -128.150053529354
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
  expect_that(lik.con(pars), equals(ll))

  ## Next, bump up alpha a little.  This really only works because
  ## theta should be basically zero I think.
  pars <- c(1, .1, 0)
  ll <- -127.649717292364
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
  expect_that(lik.con(pars), equals(ll))
})

## Here is the idea for matt;
##
## 1. Implement a prune+rescale method that won't take theta; test
## against vcv for correctness.
##
## 2. Compare what the actual calculations are (branch-for-branch) on
## a prune+rescale vs prune(ou) calculation are; that will show
## exactly what is being assumed for theta!  Do that on a two and
## three taxon tree.  Mathematica may help.
##
## 3. That will probably help add theta to the vcv and contrasts
## methods.
##
## 4. This also means that we need to implement the rescaling
## functions that arbutus does.  Try and do this efficiently and
## nicely, but we can probably just grab the functions from there.
##
## 5. That suggests that a 'bm.rescale' is possibly the correct
## interface, and that ou, method="rescale+pruning" would then call
## bm.rescale, perhaps?


## So, what happens on a three taxon tree?
str <- "((a:1,b:1):1.5,c:2.5);"
phy <- read.tree(text=str)

phy2 <- arbutus:::model.phylo.ou(phy, list(alpha=.1, sigsq=1, SE=0))




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

## These actually all disagree in their theta parameter, which is more
## than a little annoying.
expect_that(fit1$lnLik, equals(fit2$lnLik))
expect_that(coef(fit1)[1:2], equals(coef(fit2)[1:2]))

expect_that(fit1$lnLik, equals(fit3$lnLik))
expect_that(coef(fit1)[1:2], equals(coef(fit3)[1:2]))

fit4 <- no.stdout(fitContinuous(bird.orders, x, model="OU"))
expect_that(fit4$opt$lnL, equals(fit1$lnLik))
## These are quite different, but the likelihood there suggests a
## ridge...
## expect_that(c(fit4$opt$sigsq,fit4$opt$alpha),
##             equals(fit1$par[1:2], check.attributes=FALSE))
