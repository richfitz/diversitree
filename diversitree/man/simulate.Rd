\name{simulate}
\alias{trees}
\alias{tree.bd}
\alias{tree.bisse}
\alias{tree.bisseness}
\alias{tree.musse}
\alias{tree.musse.multitrait}
\alias{tree.quasse}
\alias{tree.yule}
\alias{prune}

\title{Evolve Birth-Death Trees}

\description{Evolves one or more trees under the BiSSE (Binary State
  Speciation and Extinction) birth-death model, MuSSE (Multi-State
  Speciation and Extinction) model, or a simple character independent
  birth-death model.  For BiSSE and MuSSE trees, it simultaneously
  evolves a binary character that affects speciation and/or extinction,
  and the tree itself.
}

\usage{
trees(pars, type=c("bisse", "bisseness", "bd", "musse", "quasse", "yule"),
      n=1, max.taxa=Inf, max.t=Inf, include.extinct=FALSE, ...)

tree.bisse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
           x0=NA)
tree.musse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
           x0=NA)
tree.musse.multitrait(pars, n.trait, depth, max.taxa=Inf, max.t=Inf,
                      include.extinct=FALSE, x0=NA)

tree.quasse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE, x0=NA,
            single.lineage=TRUE, verbose=FALSE)

tree.bisseness(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
               x0=NA)

tree.bd(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE)
tree.yule(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE)

prune(phy, to.drop=NULL)
}

\arguments{
  \item{pars}{Vector of parameters.  The parameters must be in the same
    order as an unconstrained likelihood function returned by
    \code{make.x}, for tree type {x}.  The MuSSE simulator automatically
    detects the appropriate number of states, given a parameter vector.}
  
  \item{type}{Type of tree to generate: May be "bisse" or "bd".}

  \item{n}{How many trees to generate?}

  \item{max.taxa}{Maximum number of taxa to include in the tree.  If
    \code{Inf}, then the tree will be evolved until \code{max.t} time
    has passed.}

  \item{max.t}{Maximum length to evolve the phylogeny over.  If
    \code{Inf} (the default), then the tree will evolve until
    \code{max.taxa} extant taxa are present.}

  \item{include.extinct}{Logical: should extinct taxa be included in
    the final phylogeny?  And should extinct trees be returned by
    \code{trees}?}
  
  \item{x0}{Initial character state at the root (state 0 or 1).  A value
    of \code{NA} will randomly choose a state from the model's
    equilibrium distribution for a BiSSE model, but a non-NA value must
    be specified for MuSSE and QuaSSE.}

  \item{n.trait, depth}{For \code{tree.musse.multitrait} only, these
    specify the number of binary traits and the style of parameters
    (with the same meaning as in \code{\link{make.musse.multitrait}}).
    The \code{pars} argument then needs to be in the same order as a
    likelihood function created by \code{make.musse.multitrait} with
    these arguments (this interface may be improved in future -- email
    me if you find this annoying).}

  \item{single.lineage}{(\code{tree.quasse} only): Start simulation with
    a single lineage?  If \code{FALSE}, then the simulation starts with
    two lineages in state \code{x0} (i.e., immediately following a
    speciation event).}
  
  \item{verbose}{(\code{tree.quasse} only): print verbose details about
    tree simuations.  This can be reassuring for really large trees.}
  
  \item{...}{Additional arguments}

  \item{phy}{A phylogeny, possibly with extinct species, produced by one
    of the tree evolving functions.}
  
  \item{to.drop}{Optional vector with the species names to drop.}
}

\details{
  The phylogeny will begin from a single lineage in state \code{x0}, but
  the final phylogeny will include only branches above the first split.

  \code{tree.bisse} may return an extinct phylogeny, and \code{trees}
  might return extinct phylogenies if \code{include.extinct} is
  \code{TRUE}.
}

\note{
  There are some logic problems around the creation of zero and one
  species trees; this will cause occasional errors when running the
  above functions.  Things will change to fix this soon.  All these
  functions may change in the near future.
}

\value{
  A \code{phylo} phylogenetic tree (ape format), or for
  \code{bisse.trees}, a \code{list} of \code{phylo} trees.

  The trees will have an element \code{tip.state} that contains the
  binary state information.
}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(3)
phy <- tree.bisse(pars, max.taxa=30, x0=0)
phy$tip.state

h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy)

## Retain extinct species:
set.seed(3)
phy2 <- tree.bisse(pars, max.taxa=30, x0=0, include.extinct=TRUE)
h2 <- history.from.sim.discrete(phy2, 0:1)
plot(h2, phy2)

#### MuSSE:
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

\keyword{models}
\author{Richard G. FitzJohn}
