## This file contains all additional functions necessary for running
## BiSSE-ness in R, as described in Magnuson-Ford, K., and
## S. Otto. 2012. Linking the investigations of character evolution
## and species diversification.  American Naturalist, XX: XX-XX.

## The BiSSE-ness model is identical to BiSSE in all aspects other
## than allowing character shifts at speciation.  This affects the
## following functions within diversitree.  The BiSSE analog of these
## functions have the same function name as those found here, omitting
## 'ness'.

## Using BiSSE-ness to estimate parameters (numbers follow other
## diversitree files)
##   1. make.bisseness is the overarching function that makes the
##      model
##
##   2. print.bisseness is the print 
##
##   3. argnames.bisseness gets and sets argument names
##
##   4. find.mle.bisseness alters the MLE search
##
##  (5. would be make.cache, but make.cache.bisse used)
##
##  (6. would be ll.bisseness, but done internally)
##
##   7. initial.conditions.bisseness combines the values from two
##      descending branches at a node, simultaneously accounting for
##      the possibility of cladogenetic change (eqn. 2 in the paper).
##
##   8. make.branches.bisseness computes calculation along a branch
##      using the differential equations specified in C file
##      'bisseness-eqs.c' (eqn. 3-6 in the paper).
##
##   9. branches.unresolved.bisseness calls the fortran code for
##      calculating the rate matrix and its exponent, to obtain the
##      probability of seeing the data listed in unresolved.
##
##  10. root.xxsseness specifies the calculations at the root, by
##      default conditioning the likelihood of the data on survival of
##      the two lineages descending from the root in a manner that
##      accounts for cladogenesis (eqn. 7; see also Nee et. al.,
##      1994).

## Using BiSSE-ness to simulate phylogenies:
##
##  11. stationary.freq.bisseness used for tree simulator to determine
##      the equilibrium root state; it can also be used to combine
##      probabilities at the root for parameter estimation using ML or
##      MCMC (but not the default, specify this explicitly).
##
##  12. tree.bisseness is the overall function used to simulate
##      phylogenies under BiSSE-ness.
##
##  13. make.tree.bisseness simulates phylogenies under BiSSE-ness.

## 1: make
make.bisseness <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                           nt.extra=10, strict=TRUE, control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend

  ## Note that bisseness uses *bisse*'s cache.
  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra,
                            strict=strict)

  if ( backend == "CVODES" )
    all.branches <- make.all.branches.C.bisseness(cache, control)
  else
    branches <- make.branches.bisseness(cache, control)

  if ( backend == "CVODES" && !is.null(cache$unresolved) )
    stop("Cannot yet use CVODES backend with unresolved clades")

  if ( !is.null(cache$unresolved) )
    warning("BiSSEness with unresolved clades has not yet been extensively tested")

  ll.bisseness <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                           root.p=NULL, intermediates=FALSE) {
    check.pars.bisseness(pars)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    if ( !is.null(cache$unresolved) )
      cache$preset <- 
        branches.unresolved.bisseness(pars, cache$unresolved)

    if ( backend == "CVODES" ) {
      ## It *should* be as simple as this, but the root treatment
      ## needs to be done differently.
      ## ll.xxsse.C(pars, all.branches,
      ##           condition.surv, root, root.p, intermediates)
      stop("Not yet written")
    } else {
      ans <- all.branches.matrix(pars, cache,
                                 initial.conditions.bisseness,
                                 branches)
      vals <- ans$init[,cache$root]
      ## This is a bit of a back, but allows us to use the
      ## root.p.xxsse function and not write our own.
      if ( root == ROOT.EQUI ) {
        root <- ROOT.GIVEN
        root.p <- stationary.freq.bisseness(pars)
        root.p <- c(root.p, 1-root.p)
      }
      root.p <- root.p.xxsse(vals, pars, root, root.p)
      loglik <- root.bisseness(vals, pars, ans$lq, condition.surv,
                               root.p)
      if ( intermediates ) {
        ans$root.p <- root.p
        attr(loglik, "intermediates") <- ans
        attr(loglik, "vals") <- vals
      }
      
      loglik
    }
  }

  class(ll.bisseness) <- c("bisseness", "function")
  ll.bisseness
}

## 2: print
print.bisseness <- function(x, ...) {
  cat("BiSSE-ness likelihood function:\n")
  print(unclass(x))
}

## 3: argnames 
argnames.bisseness <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10",
      "p0c", "p0a","p1c","p1a")
  else
    ret
}
`argnames<-.bisseness` <- function(x, value) {
  if ( length(value) != 10 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.bisseness <- function(func, x.init, method,
                           fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method,
             class.append="fit.mle.bisseness")
}

mcmc.bisseness <- mcmc.lowerzero

## 5: make.cache (make.cache.bisse used directly)

## 6: ll

## 7: initial.conditions:
## Note that we ignore both 't' and 'idx'.
initial.conditions.bisseness <- function(init, pars, t, idx) {
  lambda <- pars[1:2]
  c(init[c(1,2),1],
    init[c(3,4),1] * init[c(3,4),2] * lambda * (1 - pars[c(7,9)]) +
    ##
    ((init[c(4,4),1] * init[c(3,3),2])/2 +
     (init[c(3,3),1] * init[c(4,4),2])/2) * lambda *
    c(pars[7]*pars[8], pars[9]*pars[10]) +
    ##
    init[c(4,3),1] * init[c(4,3),2] * lambda *
    c(pars[7]*(1-pars[8]), pars[9]*(1-pars[10])))
}

## 8: branches
make.branches.bisseness <- function(cache, control) {
  neq <- 4L
  np <- 10L
  comp.idx <- as.integer(3:4)
  make.ode.branches("bisseness", "diversitree", neq, np, comp.idx,
                    control)
}

make.all.branches.C.bisseness <- function(cache, control) {
  neq <- 4L
  np <- 10L
  comp.idx <- as.integer(3:4)
  make.all.branches.C(cache, "bisseness", "diversitree",
                      neq, np, comp.idx, control)
}

## 9: branches.unresolved
branches.unresolved.bisseness <- function(pars, unresolved) {
 Nc <- unresolved$Nc
 k <- unresolved$k
 nsc <- unresolved$nsc
 t <- unresolved$len
 nt <- max(Nc) + unresolved$nt.extra

 lambda0 <- pars[1]
 lambda1 <- pars[2]
 mu0 <- pars[3]
 mu1 <- pars[4]
 q01 <- pars[5]
 q10 <- pars[6]
 p0c <- pars[7]
 p0a <- pars[8]
 p1c <- pars[9]
 p1a <- pars[10]
 base <- nucexpl(nt, lambda0, lambda1, mu0, mu1, q01, q10, p0c, p0a, p1c, p1a, t,
                 Nc, nsc, k)[,c(3,4,1,2),drop=FALSE]

 q <- rowSums(base[,3:4,drop=FALSE])
 base[,3:4] <- base[,3:4] / q

 ## Note the transpose here.
 list(target=unresolved$target,
      lq=log(q),
      base=t(base))
}

## Additional functions:
## 10. root.xxsseness
root.bisseness <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)
  is.root.both <- is.null(root.p)

  lambda <- pars[1:2]
  e.root <- vals[1:2]
  d.root <- vals[3:4]
  
  if ( condition.surv ) {
    ##  d.root <- d.root / (lambda * (1-e.root)^2) # old and incorrect
    ##  d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2) #for BiSSE
    nonextinct <-
      c(((1-pars[7])          * (1-e.root[1])^2+
         pars[7]*pars[8]      * (1-e.root[1])*(1-e.root[2]) +
         pars[7]*(1-pars[8])  * (1-e.root[2])^2),
        ((1-pars[9])          * (1-e.root[2])^2+
         pars[9]*pars[10]     * (1-e.root[1])*(1-e.root[2]) +
         pars[9]*(1-pars[10]) *(1-e.root[1])^2))
    
    d.root <- d.root / sum(root.p * lambda * nonextinct)
  }

  if ( is.root.both ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

## Additional functions for simulating trees under the BiSSE-ness model.
## 11. stationary.freq.bisseness
stationary.freq.bisseness <- function(pars) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  p0c <- pars[7]
  p0a <- pars[8]
  p1c <- pars[9]
  p1a <- pars[10]

  g   <- (lambda0 - mu0) - (lambda1 - mu1)
  eps <- (lambda0 + mu0  +  lambda1 + mu1) * 1e-14
  ## ss=state shift: i.e. change within lineages *and* at speciation.
  ss0 <- q01 + lambda0*p0c*(2*(1 - p0a) + p0a)
  ss1 <- q10 + lambda1*p1c*(2*(1 - p1a) + p1a)
  
  if ( abs(g) < eps ) {
    if ( ss0 + ss1 == 0 ) 
      0.5
    else
      ss1/(ss0 + ss1)
  } else {
    roots <- quadratic.roots(g, ss1 + ss0 - g, -ss1)
    roots <- roots[roots >= 0 & roots <= 1]
    if ( length(roots) > 1 )
      NA
    else
      roots
  }
}

## 12. tree.bisseness
tree.bisseness <- function(pars, max.taxa=Inf, max.t=Inf,
                           include.extinct=FALSE, x0=NA) {
  check.pars.bisseness(pars)
  if ( is.na(x0) )
    x0 <- as.integer(runif(1) > stationary.freq.bisseness(pars))
  else if (length(x0) != 1 || !(x0 == 0 || x0 == 1)) 
    stop("Invalid root state")

  info <- make.tree.bisseness(pars, max.taxa, max.t, x0)
  phy <- me.to.ape.bisse(info[-1, ], info$state[1])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

## 13. make.tree.bisseness
make.tree.bisseness <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                                single.lineage=TRUE)  {
  p.pars<- matrix(c(1-pars[c(7,9)],
                    pars[c(7,9)]*pars[c(8,10)],
                    pars[c(7,9)]*(1-pars[c(8,10)])), 2, 3,
                  byrow=FALSE)
  pars <- matrix(pars[1:6], 2, 3)
  
  extinct <- FALSE
  split <- FALSE
  parent <- 0
  n.i <- c(0, 0)
  r.i <- rowSums(pars)
  len <- 0
  t <- 0
  hist <- list()
  
  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0 + 1] <- 1
    start <- 0
  } else {
    stop("Nope.")
  }
  
  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    ## When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt
    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Proceed.  What state does an event happen to?
    state <- as.integer(runif(1) > r.n[1]/r.tot)
    state.i <- state + 1

    ## Pick a lineage for that state:
    j <- sample(n.i[state.i], 1)
    lineage <- lineages[states[lineages] == state][j]

    ## Pick an event: 1= speciation, 2= extinction, 3= state change
    type <- sample(3, 1, FALSE, pars[state.i, ])

    if ( type == 1 ) {
      ## Speciating:      
      if (n.taxa == max.taxa) 
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      
      if ( p.pars[state.i,1]==1 ) {
        states[new.i] <- state
        n.i[state.i] <- n.i[state.i] + 1
      } else {
        ## The daughter states be inherited from the parent state? 
        ## 1= complete inheritance (i.e. traditional bisse), 
        ## 2= partial inheritance, 
        ## 3= no inheritance.
        inherit.type<- sample(3, 1, FALSE, p.pars[state.i,])

        ## RGF: I suspect that moving the statements from this switch
        ## into the three-way if below will be more efficient.
        states[new.i] <- switch(inherit.type,
                                state * c(1,1),
                                {new1<-sample(0:1, 1); c(new1, 1-new1)},
                                rep(1-state, 2))
        
        if ( inherit.type == 1 )
          n.i[state.i] <- n.i[state.i] + 1
        else if ( inherit.type == 2 )
          n.i[(1-state)+1] <- n.i[(1-state)+1] + 1
        else if ( inherit.type == 3 ) {
          n.i[(1-state)+1] <- n.i[(1-state)+1] + 2
          n.i[state.i] <- n.i[state.i] - 1
        }
      }
      
      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      n.taxa <- n.taxa + 1
      lineages <- which(!split & !extinct)
    } else if ( type == 2 ) {
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.i[state.i] <- n.i[state.i] - 1
      n.taxa <- n.taxa - 1
    } else {
      n.i <- n.i + if (state == 0) c(-1, 1) else c(1, -1)
      states[lineage] <- 1 - state
      hist[[length(hist) + 1]] <- c(lineage, t, state, 1 - state)
    }
  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent, 
                     start=start, state=states, extinct=extinct,
                     split=split)
  hist <- as.data.frame(do.call(rbind, hist))
  if (nrow(hist) == 0)
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0

  attr(info, "t") <- t
  attr(info, "hist") <- hist
  info
}

check.pars.bisseness <- function(pars) {
  if ( length(pars) != 10 )
    stop("Invalid parameter length (expected 10)")
  if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
  ## Already checked "< 0" above.
  if ( any(pars[7:10] > 1) )
    stop("Probability parameters must lie between 0 and 1.")
  TRUE
}

