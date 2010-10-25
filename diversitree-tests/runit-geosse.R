## Test case code, based on Emma's tests.
test.geosse <- function() {
  tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")

  states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
  names(states) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")

  starting.point.geosse(tree)
                                        #        sA        sB       sAB        xA        xB        dA        dB 
                                        # 1.5788510 1.5788510 1.5788510 0.0000000 0.0000000 0.3157702 0.3157702 

  pars <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)
  lnL.7par <- make.geosse(tree, states)
  names(pars) <- argnames(lnL.7par)

  checkEquals(lnL.7par(pars), -23.71574, tolerance=1e-7)
  checkEquals(lnL.7par(pars, condition.surv=TRUE), -24.10259,
              tolerance=1.06e-7)

  lnL <- constrain(lnL.7par, sA ~ sB)
  checkEquals(lnL(pars[2:7]), -24.34399, tolerance=1e-7)
  checkEquals(lnL(pars[2:7], condition.surv=TRUE), -24.6769,
              tolerance=1e-7)

  lnL <- constrain(lnL.7par, sAB ~ 0)
  checkEquals(lnL(pars[-3]), -23.56325, tolerance=1e-7)
  checkEquals(lnL(pars[-3], condition.surv=TRUE), -23.90416,
              tolerance=2e-7)

  lnL <- constrain(lnL.7par, dB ~ dA)
  checkEquals(lnL(pars[-6]), -23.2386, tolerance=1e-7)

  lnL <- constrain(lnL.7par, dB ~ dA, sAB ~ 0)
  checkEquals(lnL(pars[-c(3,7)]), -23.763)

  fit <- find.mle(lnL.7par, pars)
  checkEquals(fit$lnLik, -18.77671, tolerance=2e-7)

  fit <- find.mle(constrain(lnL.7par, dA ~ dB), pars[-6])
  checkEquals(fit$lnLik, -18.77672, tolerance=2e-7)

  fit <- find.mle(constrain(lnL.7par, dA ~ dB, sAB ~ 0), pars[-c(3,7)])
  checkEquals(fit$lnLik, -18.77672, tolerance=2e-7)

  checkEquals(lnL.7par(pars, root.p=c(0.9, 0.1, 0), root=ROOT.GIVEN),
              -23.56753, tolerance=1e-7)
  checkEquals(lnL.7par(pars, root=ROOT.EQUI), -24.32152, tolerance=1e-7)
  checkEquals(lnL.7par(pars, root=ROOT.FLAT), -24.25945, tolerance=2e-7)

  set.seed(1)
  samples <- mcmc(lnL.7par, pars, nsteps=3, lower=0, upper=3, w=3/5,
                  prior=make.prior.exponential(1)) 
  checkEquals(samples$p[3], -27.29315, tolerance=2e-7)

  set.seed(1)
  samples <- mcmc(constrain(lnL.7par, dA ~ dB, sAB ~ 0), pars[-c(3,7)],
                  nsteps=3, lower=0, upper=3, w=3/5,
                  prior=make.prior.exponential(1))
  checkEquals(samples$p[3], -25.94858)
}
