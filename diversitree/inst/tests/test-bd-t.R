## TODO: Could do with more extensive tests here.
library(diversitree)
library(testthat)

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("Birth-death (time-dependent)")

set.seed(1)
pars <- c(.1, .03)
phy <- trees(pars, "bd", max.taxa=25)[[1]]

## Next, make three different likelihood functions: a "normal" one that
## uses the direct birth-death calculation, an "ode" based one (that
## uses numerical integration to compute the likelihood, and is
## therefore not exact), and one that is time-varying, but that the
## time-dependent functions are constant.t().
lik.nee <- make.bd(phy)
lik.ode <- make.bd(phy, control=list(method="ode"))
lik.t <- make.bd.t(phy, c("constant.t", "constant.t"))

expect_that(lik.nee(pars), equals(-22.50266563415211))

## ODE-based likelihood calculations are correct to about 1e-6.
expect_that(lik.ode(pars), equals(lik.nee(pars)))

## The ODE calculation agrees exactly with the time-varying (but
## constant) calculation.
expect_that(lik.ode(pars), is_identical_to(lik.t(pars)))

## Next, make a real case, where speciation is a linear function of
## time.
lik.t2 <- make.bd.t(phy, c("linear.t", "constant.t"))

## Confirm that this agrees with the previous calculations when the
## slope is zero
pars2 <- c(pars[1], 0, pars[2])
expect_that(lik.t2(pars2),  is_identical_to(lik.t(pars)))

set.seed(1)
pars3 <- pars2 + runif(length(pars2), 0, .2)
expect_that(lik.t2(pars3), equals(-105.437123043537))
