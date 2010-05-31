\name{tree.musse}
\alias{tree.musse}
\alias{trees}

\title{Evolve Birth-Death Trees Under MuSSE}

\description{Evolves a tree under the MuSSE (Multi State Speciation and
  Extinction) birth-death model, simultaneously evolving a discrete
  character that affects speciation and/or extinction, and the tree
  itself.
}

\usage{
make.tree.musse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
                x0=NA, ...)
}

\arguments{
  \item{pars}{Vector of parameters; see Details for order.}
  \item{max.taxa}{Maximum number of taxa to include in the tree.  If
    \code{Inf}, then the tree will be evolved until \code{max.t} time
    has passed.}
  \item{max.t}{Maximum length to evolve the phylogeny over.  If
    \code{Inf} (the default), then the tree will evolve until
    \code{max.taxa} extant taxa are present.}
  \item{include.extinct}{Logical: should extinct taxa be included in
    the final phylogeny?  And should extinct trees be returned by
    \code{trees}?}
  \item{x0}{Initial character state at the root (state 1 to \code{k},
    \code{k} being the highest state code).}
  \item{...}{Additional arguments, ignored.}
}


\details{
  The phylogeny will begin from a single lineage in state \code{x0}, but
  the final phylogeny will include only branches above the first split.

  The parameters are in the same order as \code{make.musse}: For
  \code{k} states, the first \code{k} parameters are speciation rates,
  the next \code{k} parameters are speciation rates, and the remaining
  \code{k(k-1)} parameters are the Q matrix in the order
  \code{q12, ...,q1k, q21, q23, ..., q2k, ..., qk(k-1)}.  For two
  states, this is the same order as for BiSSE.

  Note that multiple states can be many traits, each of which is binary
  or multi-state.  All that matters is that the number of possible
  state/trait combinations is discrete.
}


\value{
  A \code{phylo} phylogenetic tree (ape format).

  The trees will have an element \code{tip.state} that contains the
  binary state information.  The function
  \code{\link{history.from.sim.discrete}} can be used to examine fossil
  character states (see Examples).
}

\examples{
## Two states
pars <- c(.1, .2, .03, .04, 0.05, 0.1)
set.seed(2)
phy <- tree.musse(pars, 20, x0=1)

h <- history.from.sim.discrete(phy, 1:2)
plot(h, phy)

## A 3-state example where movement is only allowed between neighbouring
## states (1 <-> 2 <-> 3), and where speciation and extinction rates
## increase moving from 1 -> 2 -> 3:
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32

set.seed(2)
phy <- tree.musse(pars, 30, x0=1, include.extinct=TRUE)

h <- history.from.sim.discrete(phy, 1:3)
plot(h, phy, cex=.7)

## And with extinct taxa pruned:
phy2 <- prune(phy)
h2 <- history.from.sim.discrete(phy2, 1:3)
plot(h2, phy2, cex=.7)

## This can all be done in one step (and is by default):
set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
h <- history.from.sim.discrete(phy, 1:3)
plot(h, phy, cex=.7)
}