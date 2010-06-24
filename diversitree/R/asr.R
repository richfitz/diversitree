## Core ancestral state reconstruction (ASR) code.  I am implementing
## three different types of things:
##   (1) asr.marginal
##   (2) asr.joint
##   (3) asr.stoch
## There will perhaps be an overarching "asr" function some day but
## not right now.  All ASR functions take as a first argument a
## likelihood function, and dispatch based on the class of this.
## Methods currently implemented are mk2/mkn, but bisse and musse
## methods can be found in the unreleased package diversitree.unrel

## Core generics
asr.marginal <- function(lik, pars, nodes=NULL, ...) {
  UseMethod("asr.marginal")
}

asr.joint <- function(lik, pars, n=1, simplify=TRUE, ...) {
  UseMethod("asr.joint")
}

asr.stoch <- function(lik, pars, n=1, ...) {
  UseMethod("asr.stoch")
}

## Method for constrained function dispatch on the next class down: a
## constrained mkn model has class:
##   c("constrained", "mkn")
## so asr.xxx dispatches to asr.xxx.constrained -> asr.xxx.mkn
## with the appropriate parameters filled in by asr.xxx.constrained.
## This trick works for all three functions:
asr.marginal.constrained <- function(lik, pars, nodes=NULL, ...) {
  pars <- lik(pars, pars.only=TRUE)
  lik <- attr(lik, "func")
  NextMethod("asr.marginal")
}

asr.joint.constrained <- function(lik, pars, n=1, simplify=TRUE, ...) {
  pars <- lik(pars, pars.only=TRUE)
  lik <- attr(lik, "func")
  NextMethod("asr.joint")
}

asr.stoch.constrained <- function(lik, pars, n=1, ...) {
  pars <- lik(pars, pars.only=TRUE)
  lik <- attr(lik, "func")
  NextMethod("asr.stoch")
}

## Next, the utility functions for the different types of models
## This is to asr.marginal what all.branches is for the core models.
## Unfortunately, it is a little unclear what one should do at the
## root, so this is (for BD-based models) *not* conditioning on
## survival and for all models is using ROOT.OBS to combine Ds at the
## root.  I will probably have to allow for a "root" function to be
## used here to get around this this though.
##
## The argument 'res' is the result of running all.branches
##
## TODO: There is a substantial optimisation for doing multiple
## samples at once; at each node, we just have to interate over the
## possible states at that node, as these will match for
##   di * pij[nd,i,]
## And then sample the appropriate number of points.
do.asr.marginal <- function(pars, cache, res, nodes, states.idx,
                            initial.conditions,
                            branches, root,
                            ...) {
  ## Store these for easier calculation.
  children <- cache$children
  parent <- cache$parent
  len <- cache$len
  depth <- cache$depth
  root.idx <- cache$root
  anc <- cache$anc

  if ( is.null(nodes) )
    nodes <- root.idx:max(cache$order)
  else
    nodes <- nodes + cache$n.tip

  f <- function(nd) {
    ## Include current node but omit root:
    anc.nd <- c(nd, anc[[nd]])
    anc.nd <- anc.nd[-length(anc.nd)]
    p <- rep(NA, length(states.idx))
    
    for ( st in seq_along(states.idx) ) {
      lq <- res$lq
      branch.init <- res$init
      branch.base <- res$base
      branch.init[[nd]][states.idx[-st]] <- 0
      y.in <- branch.init[[nd]]
      ## j <- nd # Needed for when nd == root.idx - TODO:poss no longer?

      for ( i in anc.nd ) {
        ans <- branches(y.in, len[i], pars, depth[i])
        lq[i] <- ans[1]
        branch.base[[i]] <- ans[-1]
        j <- parent[i]
        y.in <- initial.conditions(branch.base[children[j,]], pars,
                                   depth[j], j == root.idx)
        branch.init[[j]] <- y.in
      }

      ans <- root(pars, branch.init[[root.idx]], lq)
      
      if ( is.na(ans) )
        p[st] <- -Inf # explots R's exp(-Inf) == 0
      else
        p[st] <- ans
    }

    pp <- exp(p - max(p))
    pp/sum(pp)
  }

  matrix(unlist(lapply(nodes, f)), ncol=length(nodes))
}

## Utility function for drawing one or more samples from the joint
## distribution.
## TODO: Allow root.state to be either a percentage or a function of
## the parameters?  -- no need for the latter, as the parameters do
## not change in this.  root.p[i] gives the probability that the root is
## in state i; this can probably replace root.state, as the two just
## fight each other.

## li is len * k matrix; for a node n, li[n,] comes in the order
##   Pr(D_n|1), Pr(D_n|2), ..., Pr(D_n|k)
## Pr(D|i) is the conditional probability of the data conditional on a
## node being in the state 'i'.  It .

## pij is a len * (k*k) column matrix; the row n comes in the order
##   p11, p21, ..., pk1, p12, ..., pkk
## so that
##   matrix(pij[nd,], k, k)
## is a matrix with where m[i,j] is the probability of moving from
## state i to state j.  This means that
##   pij2 <- array(pij, c(nrow(pij), k, k))
## gives an array where
##   pij2[nd, i, j]
## is the probability of an i->j transition along the branch leading
## to nd.  Note that the sums of all rows (pij2[n,i,] for all n, i)
## equals 1, as a branch ends at some state with probability 1:
##   all(abs(apply(pij2, 1:2, sum)[-root,] - 1) < 1e-8)
do.asr.joint <- function(n, cache, li, pij, root.p, simplify=TRUE,
                         ...) {
  parent <- cache$parent
  len <- length(cache$len)
  root <- cache$root
  tips <- seq_len(cache$n.tip)
  k <- length(li[[1]])
  pij2 <- array(pij, c(len, k, k))

  anc.states <- matrix(NA, n, len)
  anc.states[,root] <- sample(k, n, TRUE, root.p)
  for ( i in rev(cache$order)[-1] ) {
    parent.state <- anc.states[,parent[i]]
    for ( j in seq_len(k) ) {
      idx <- parent.state == j
      nj <- sum(idx)
      if ( nj > 0 ) {
        p <- li[[i]] * pij2[i,j,] # di * pij
        anc.states[idx,i] <- sample(k, nj, TRUE, p)
      }
    }
  }
  
  anc.states[,-tips,drop=simplify]
}

## TODO: explain this one...
do.asr.jointmean <- function(cache, li, pij, root.p, ...) {
  parent <- cache$parent
  len <- length(cache$len)
  root <- cache$root
  tips <- seq_len(cache$n.tip)
  k <- length(li[[1]])
  pij2 <- array(pij, c(len, k, k))

  anc.states <- matrix(NA, k, len)  
  anc.states[,root] <- root.p

  for ( i in rev(cache$order)[-1] ) {
    pp <- anc.states[,parent[i]]
    tmp <- t((li[[i]] * t(pij2[i,,])))
    anc.states[,i] <- (c(pp) / sum(pp)) %*% (tmp/rowSums(tmp))
  }

  anc.states[,-tips]
}
