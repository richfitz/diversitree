\name{make.musse}
\alias{make.musse}
\alias{starting.point.musse}

\title{MuSSE: Multi-State Speciation and Extinction}

\description{Prepare to run MuSSE (Multi-State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.

  MuSSE is agnostic as to whether multiple states or multiple traits are
  modelled (following Pagel 1994).  Instead, a transition rate matrix
  amongst possible trait/state combinations is constructed and the
  analysis is conducted on this.
}

\usage{
make.musse(tree, states, k, sampling.f=NULL, strict=TRUE)
starting.point.bisse(tree, q.div=5, yule=FALSE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape}
    \dQuote{phylo} format.}

  \item{states}{A vector of character states,
    each of which must be 0 or 1.  This vector must have names that
    correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips corresponding to unresolved clades,
    the state should be \code{NA}.}

  \item{k}{The number of states.}
  
  \item{sampling.f}{Vector of length \code{k} where \code{sampling.f[i]}
    is the proportion of species in state \code{i} that are present in
    the phylogeny.  A value of \code{c(0.5, 0.75, 1)} means that half of
    species in state 1, three quarters of species in state 2, and all
    species in state 3 are included in the phylogeny.  By default all
    species are assumed to be known}

  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are integers on \code{1:k}.  If \code{strict} is
    \code{TRUE} (the default), then the additional check is made that
    \emph{every} state is present.  The likelihood models tend to be
    poorly behaved where states are missing, but there are cases
    (missing intermediate states for meristic characters) where allowing
    such models may be useful.}

  \item{q.div}{Ratio of diversification rate to character change rate.
    Eventually this will be changed to allow for Mk2 to be used for
    estimating q parameters.}
  
  \item{yule}{Logical: should starting parameters be Yule estimates
    rather than birth-death estimates?}
}

\examples{
## 1: BiSSE equivalence
pars <- c(.1, .2, .03, .04, 0.05, 0.1)
set.seed(2)
phy <- tree.musse(pars, 20, x0=1)

## Show that the likelihood functions give the same answers.  Ignore the
## warning when creating the MuSSE function.
lik.b <- make.bisse(phy, phy$tip.state-1)
lik.m <- make.musse(phy, phy$tip.state, 2)
all.equal(lik.b(pars), lik.m(pars), tolerance=1e-7)

## Notice that default argument names are different between BiSSE and
## MuSSE, but that the order is the same.
argnames(lik.b) # BiSSE: 0/1
argnames(lik.m) # MuSSE: 1/2

## 2: A 3-state example where movement is only allowed between
## neighbouring states (1 <-> 2 <-> 3), and where speciation and
## extinction rates increase moving from 1 -> 2 -> 3:
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
lik <- make.musse(phy, phy$tip.state, 3)

## Quite a few parameters now!  This will slow the ML search down badly.
argnames(lik)

## Start with a fully constrained model
constraints <- list(lambda2 ~ lambda1,



}

\keyword{models}