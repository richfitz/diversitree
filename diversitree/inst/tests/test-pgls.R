library(diversitree)
library(testthat)
library(caper)
library(nlme)

context("PGLS")

## Simulated tree and traits:
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
ty <- sim.character(phy, 1) + 3
ta <- sim.character(phy, 1)
tb <- sim.character(phy, 1)

data <- data.frame(a=ta, b=tb, y=ty, row.names=names(ty))

## Fit the model using standard approaches:
## 1. caper::pgls
cdata <- comparative.data(phy, cbind(data, species=phy$tip.label),
                          "species")
fit.caper.ya  <- pgls(y ~ a, cdata)
fit.caper.yab <- pgls(y ~ a + b, cdata)

## 2. nlme::gls
fit.gls.ya  <- gls(y ~ a,     data, corBrownian(1, phy), method="ML")
fit.gls.yab <- gls(y ~ a + b, data, corBrownian(1, phy), method="ML")

## This is a direct simplification of the code from Freckleton (2012);
## skipping missing data treatment, and returning the ML diffusion
## estimate.  It is this last part that we need here.
pgls.contrasts <- function(formula, data, phylo) {
  cdata <- data.frame(apply(data, 2, pic, phylo))
  fm <- lm(formula, cdata, y=TRUE)
  vars <- pic(data[[1]], phylo, var=TRUE)[,2]
  V <- diversitree:::pgls.root.var.bm(phylo)
  u <- residuals(fm) 
  n <- length(u)
  sigma2 <- sum(u^2) / (n + 1)
  logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars ))
                   +  log(V) + sum(u^2) / sigma2)
  list(model=fm, logLik=logLik, sigma2=sigma2)
}

fit.contrasts.ya  <- pgls.contrasts(y ~ a     - 1, data, phy)
fit.contrasts.yab <- pgls.contrasts(y ~ a + b - 1, data, phy)

test_that("Reference implementations agree", {
  expect_that(fit.contrasts.ya$logLik,
              equals(as.numeric(logLik(fit.gls.ya))))
  expect_that(fit.contrasts.yab$logLik,
              equals(as.numeric(logLik(fit.gls.yab))))

  expect_that(coef(fit.contrasts.ya$model),
              equals(coef(fit.caper.ya)[-1]))
  expect_that(coef(fit.contrasts.yab$model),
              equals(coef(fit.caper.yab)[-1]))

  expect_that(coef(fit.gls.ya),  equals(coef(fit.caper.ya)))
  expect_that(coef(fit.gls.yab), equals(coef(fit.caper.yab)))
})

## From this, extract the full set of ML parameters (this requires
## interrogating both the caper and contrasts fits)
p.ya <- unname(c(coef(fit.caper.ya), fit.contrasts.ya$sigma2))
p.yab <- unname(c(coef(fit.caper.yab), fit.contrasts.yab$sigma2))
l.ya <- fit.contrasts.ya$logLik
l.yab <- fit.contrasts.yab$logLik

## Arbitrary offset:
p2.ya  <- p.ya  + runif(length(p.ya),  0, 0.1)
p2.yab <- p.yab + runif(length(p.yab), 0, 0.1)

## 4. VCV equation
lik.vcv.ya  <- make.pgls(phy, y ~ a,     data, control=list(method="vcv"))
lik.vcv.yab <- make.pgls(phy, y ~ a + b, data, control=list(method="vcv"))

test_that("VCV Likelihoods agree at ML points", {
  expect_that(lik.vcv.ya(p.ya),   equals(l.ya))
  expect_that(lik.vcv.yab(p.yab), equals(l.yab))
})

## 5. Contrasts
lik.con.ya  <- make.pgls(phy, y ~ a,     data, control=list(method="contrasts"))
lik.con.yab <- make.pgls(phy, y ~ a + b, data, control=list(method="contrasts"))

test_that("Contrasts Likelihoods agree at ML points", {
  expect_that(lik.con.ya(p.ya),   equals(l.ya))
  expect_that(lik.con.yab(p.yab), equals(l.yab))
})

test_that("Calculations agree at differing point", {
  expect_that(lik.con.ya(p2.ya),   equals(lik.vcv.ya(p2.ya)))
  expect_that(lik.con.yab(p2.yab), equals(lik.vcv.yab(p2.yab)))
})
