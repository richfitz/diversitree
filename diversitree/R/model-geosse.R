## GeoSSE model, by Emma Goldberg <eeg@uic.edu>

## Models should provide:
##   1. make
##   2. print
##   3. argnames / argnames<-
##   4. find.mle
## Generally, make will require:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches

## 1: make
make.geosse <- function(tree, states, sampling.f=NULL, strict=TRUE,
                        control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend

  unresolved <- NULL
  cache <- make.cache.geosse(tree, states, sampling.f, strict)

  if ( backend == "CVODES" )
    all.branches <- make.all.branches.C.geosse(cache, control)
  else
    branches <- make.branches.geosse(cache, control)

  initial.conditions <- initial.conditions.geosse

  ll.geosse <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                        root.p=NULL, intermediates=FALSE) {
    check.pars.geosse(pars)
    if ( !is.null(root.p) && root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    ## hack (from bisseness) to avoid needing root.p.geosse()
    if ( root == ROOT.EQUI ) {
      root <- ROOT.GIVEN
      root.p <- stationary.freq.geosse(pars)
    }

    if ( backend == "CVODES" ) {
      ans <- all.branches(pars, intermediates)
      lq <- ans[[1]]
      vals <- ans[[2]]
      ans <- ans$intermediates
    } else {
      ans <- all.branches.matrix(pars, cache,
                                 initial.conditions, branches)
      lq <- ans$lq
      vals <- ans$init[,cache$root]
    }
    root.p <- root.p.xxsse(vals, pars, root, root.p)
    loglik <- root.geosse(vals, pars, lq, condition.surv, root.p)

    if ( intermediates ) {
      ans$root.p <- root.p
      attr(loglik, "intermediates") <- ans
      attr(loglik, "vals") <- vals
    }
    loglik
  }

  class(ll.geosse) <- c("geosse", "function")
  ll.geosse
}

## 2: print
print.geosse <- function(x, ...) {
  cat("GeoSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.geosse <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
  else
    ret
}
`argnames<-.geosse` <- function(x, value) {
  if ( length(value) != 7 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.geosse <- function(func, x.init, method,
                            fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.geosse")
}

mcmc.geosse <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
## almost identical to make.cache.bisse(), but uses three states
make.cache.geosse <- function(tree, states, sampling.f=NULL, 
                              strict=TRUE) {
  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=strict, strict.vals=0:2)

  states <- check.integer(states)
  sampling.f <- check.sampling.f(sampling.f, 3)
  
  cache <- make.cache(tree)
  cache$ny <- 6L
  cache$tip.state  <- states
  cache$sampling.f <- sampling.f
  cache$y <- initial.tip.geosse(cache)

  cache
}

## By the time this hits, unresolved clades and any other non-standard
## tips have been removed.  We have an index "tips" (equal to 1:n.tip
## for plain bisse/geosse) that is the "index" (in phy$edge numbering) of the
## tips, and a state vector cache$tip.state, both of the same length.
## The length of the terminal branches is cache$len[cache$tips].
##
## Initial conditions at the tips are given by their tip states:
## There are four types of initial condition in geosse:
##             E0   E1   E2    D0     D1     D2
##   state0: c(1-f_0, 1-f_1, 1-f_2, f_0, 0,   0, )
##   state1: c(1-f_0, 1-f_1, 1-f_2, 0,   f_1, 0, )
##   state2: c(1-f_0, 1-f_1, 1-f_2, 0,   0,   f_2)
##   state?: c(1-f_0, 1-f_1, 1-f_2, f_0, f_1, f_2)
initial.tip.geosse <- function(cache) {
  f <- cache$sampling.f
  y <- list(c(1-f, f[1], 0, 0),
            c(1-f, 0, f[2], 0),
            c(1-f, 0, 0, f[3]),
            c(1-f, f))
  y.i <- cache$tip.state + 1
  y.i[is.na(y.i)] <- 4

  tips <- cache$tips

  ## This would return a data structure appropriate for the more basic
  ## tip treatment:
  ##   dt.tips.ordered(y[y.i], tips, cache$len[tips])
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}

## 6: ll

## 7: initial.conditions:
initial.conditions.geosse <- function(init, pars, t, idx) {
  ## E.0, E.1, E.2
  e <- init[c(1,2,3),1]

  ## D.1, D.2  (Eq. 6bc)
  d12 <- init[c(5,6),1] * init[c(5,6),2] * pars[c(1,2)]

  ## D.0 (Eq. 6a)
  d0 <- 0.5 * sum(init[c(4,5),1] * init[c(5,4),2] * pars[1] + 
                  init[c(4,6),1] * init[c(6,4),2] * pars[2] +
                  init[c(5,6),1] * init[c(6,5),2] * pars[3])
  d <- c(d0, d12)

  c(e, d)
}

## 8: branches
make.branches.geosse <- function(cache, control) {
  neq <- 6L
  np <- 7L
  comp.idx <- as.integer(4:6)
  make.ode.branches("geosse", "diversitree", neq, np, comp.idx,
                    control)
}

make.all.branches.C.geosse <- function(cache, control) {
  neq <- 6L
  np <- 7L
  comp.idx <- as.integer(4:6)
  make.all.branches.C(cache, "geosse", "diversitree",
                      neq, np, comp.idx, control)
}

## Additional functions
stationary.freq.geosse <- function(pars) {
  sA  <- pars[1]
  sB  <- pars[2]
  sAB <- pars[3]
  xA  <- pars[4]
  xB  <- pars[5]
  dA  <- pars[6]
  dB  <- pars[7]

  A <- matrix(c(
               -xA - xB - sAB,   dA,             dB,
                sA + xB + sAB,   sA - dA - xA,   0,
                sB + xA + sAB,   0,              sB - dB - xB
               ), nrow=3, byrow=TRUE)

  ## continuous time, so the dominant eigenvalue is the largest one
  ## return its eigenvector, normalized  
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  evA$vectors[,i] / sum(evA$vectors[,i])
}

## Replacement function from Emma, 4 Nov 2010.
## from Magallon & Sanderson (2001), rather than a bd fit
starting.point.geosse <- function(tree, eps=0.5) {
  if (eps == 0) {
    s <- (log(Ntip(tree)) - log(2)) / max(branching.times(tree))
    x <- 0
    d <- s/10
  } else {
    n <- Ntip(tree)
    r <- ( log( (n/2) * (1 - eps*eps) + 2*eps + (1 - eps)/2 *
               sqrt( n * (n*eps*eps - 8*eps + 2*n*eps + n))) - log(2)
          ) / max(branching.times(tree))
    s <- r / (1 - eps)
    x <- s * eps
    d <- x
  }
  p <- c(s, s, s, x, x, d, d)
  names(p) <- argnames.geosse(NULL)
  p
}

## modified from diversitree-branches.R: root.xxsse()
##   the only difference is lambda
root.geosse <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)

  k <- length(vals) / 2
  i <- seq_len(k)

  e.root <- vals[i]
  d.root <- vals[-i]

  if ( condition.surv )
  {
    ## AB species are subject to all three speciation rates
    lambda <- c(sum(pars[1:3]), pars[1:2])
    ## d.root <- d.root / (lambda * (1-e.root)^2) # old
    d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2)
  }

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

check.pars.geosse <- function(pars) {
  if ( length(pars) != 7 )
    stop("Invalid parameter length (expected 7)")
  if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
  TRUE
}
