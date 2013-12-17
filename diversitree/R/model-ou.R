## TODO: It should be fairly straightforward to rewrite a pruning
## calculator that uses the rescaling approach; we are just dealing
## with a variance of eq:3 in ou-3.tex.  Super easy.  Just sort out
## how the indexing works and all should come together.

## That means that we probably need an additional argument to the
## function: "free.optimum=TRUE/FALSE"; that should work.  Check that
## the optimum *can* be freed (only with pruning) and work from there.

## My simple-minded OU calculator, direct from the SDE:

## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches

## 1: make
make.ou <- function(tree, states, states.sd=0, control=list()) {
  control <- check.control.continuous(control)
  cache <- make.cache.ou(tree, states, states.sd, control)

  if (control$method == "vcv") {
    all.branches <- make.all.branches.ou.vcv(cache, control)
    rootfunc <- rootfunc.bm.vcv
  } else if (control$method == "pruning") {
    all.branches <- make.all.branches.ou.pruning(cache, control)
    rootfunc <- rootfunc.bm.pruning
  } else if (control$method == "contrasts") {
    all.branches <- make.all.branches.ou.contrasts(cache, control)
    rootfunc <- rootfunc.bm.contrasts
  } else {
    stop("Unknown method", method)
  }

  ll <- function(pars, root=ROOT.MAX, root.x=NULL,
                 intermediates=FALSE) {
    check.pars.ou(pars)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, root, root.x, intermediates)
  }
  class(ll) <- c("ou", "dtlik", "function")
  ll
}

## 2: info
make.info.ou <- function(phy) {
  list(name="ou",
       name.pretty="Ornstein-Uhlenbeck",
       ## Parameters:
       np=3L,
       argnames=default.argnames.ou(),
       ## Variables:
       ny=3L,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "I really don't know"))
}
default.argnames.ou <- function() c("s2", "alpha", "theta")

## 3: make.cache
make.cache.ou <- function(tree, states, states.sd, control) {
  cache <- make.cache.bm(tree, states, states.sd, control)
  cache$info <- make.info.ou(tree)
  cache
}

###########################################################################
## Additional functions

## Checking
check.pars.ou <- function(pars) {
  if ( length(pars) != 3 )
    stop("Incorrect parameter length")
  check.nonnegative(pars[1:2])
  TRUE
}
