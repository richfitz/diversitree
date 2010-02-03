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

asr.joint <- function(lik, pars, n=1, root.state=NA, simplify=TRUE,
                      ...) {
  UseMethod("asr.joint")
}

asr.stoch <- function(lik, pars, n=1, root.state=NA, ...) {
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

asr.joint.constrained <- function(lik, pars, n=1, root.state=NA, ...) {
  pars <- lik(pars, pars.only=TRUE)
  lik <- attr(lik, "func")
  NextMethod("asr.joint")
}

asr.stoch.constrained <- function(lik, pars, n=1, root.state=NA, ...) {
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
do.asr.marginal <- function(pars, cache, res, nodes, states.idx,
                            initial.conditions,
                            branches, branches.unresolved, ...) {
  ## Store these for easier calculation.
  children <- cache$children
  parent <- cache$parent
  len <- cache$len
  depth <- cache$depth
  root <- cache$root
  anc <- cache$anc

  if ( is.null(nodes) )
    nodes <- root:max(cache$order)
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
      branch.init[nd,states.idx[-st]] <- 0
      y.in <- branch.init[nd,]
      j <- nd # Needed for when nd == root

      for ( i in anc.nd ) {
        ans <- branches(y.in, len[i], pars, depth[i])
        lq[i] <- ans[1]
        branch.base[i,] <- ans[-1]
        j <- parent[i]
        y.in <- initial.conditions(branch.base[children[j,],], pars,
                                   depth[j], j == root)
        branch.init[j,] <- y.in
      }

      ## This bit is the root calculation:
      ## TODO: function that takes
      ##   branch.init[cache$root,]
      ## and the parameters, converting it into a likelihood?
      d.root <- branch.init[cache$root,states.idx]
      if ( sum(d.root) > 0 ) {
        p.root <- d.root/sum(d.root) # ROOT.OBS
        p[st] <- log(sum(p.root * d.root)) + sum(lq)
      } else {
        p[st] <- -Inf # explots R's exp(-Inf) == 0
      }
    }

    pp <- exp(p - max(p))
    pp/sum(pp)
  }

  matrix(unlist(lapply(nodes, f)), ncol=length(nodes))
}

## Utility function for drawing one or more samples from the joint
## distribution.
do.asr.joint <- function(pars, n, root.state, cache, li, pij,
                         node.labels=NULL, simplify=TRUE, ...) {
  parent <- cache$parent
  len <- length(cache$len)
  root <- cache$root
  tips <- seq_len(cache$n.tip)
  
  f <- function() {
    k <- ncol(li)
    idx <- matrix(seq_len(k*k), k, k)
    anc.states <- rep(as.numeric(NA), len)

    if ( is.na(root.state) )
      root.state <- sample(k, 1, FALSE, li[root,])
    anc.states[root] <- root.state
    
    for ( i in rev(cache$order)[-1] ) {
      parent.state <- anc.states[parent[i]]
      p <- li[i,] * pij[i,idx[,parent.state]] # di * pij
      anc.states[i] <- sample(k, 1, FALSE, p)
    }

    structure(anc.states[-tips], names=node.labels)
  }
  
  if ( n == 1 && simplify )
    x <- f()
  else {
    x <- replicate(n, f(), simplify)
    if ( simplify )
      x <- t(x)
  }

  x
}
