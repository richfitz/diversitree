\name{simulate}
\alias{trees}
\alias{tree.bisse}
\alias{tree.bd}

\title{Evolve Birth-Death Trees}

\description{Evolves one or more trees under the BiSSE (Binary State
  Specialtion and Extinction) birth-death model or a simple character
  independent birth-death model.  For BiSSE trees, it simultaneously
  evolves a binary character that affects speciation and/or extinction
  and lineages.
}

\usage{
trees(pars, type=c("bisse", "bd"), n=1, max.taxa=Inf,
      max.t=Inf, include.extinct=FALSE, ...)
tree.bisse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
           x0=NA, ...)
tree.bd(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE, ...)
}

\arguments{
  \item{pars}{Vector of parameters, in the order \code{lambda0},
    \code{lambda1}, \code{mu0}, \code{mu1}, \code{q01}, \code{q10} for a
    BiSSE tree, and in the order \code{lambda}, \code{mu} for a
    birthdeath tree.  Note that this is the same order as the function
    returned by \code{\link{make.bisse}}.}
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
    equilibrium distribution.}
  \item{...}{Additional arguments}
}

\details{
  The phylogeny will begin from a single lineage in state \code{x0}, but
  the final phylogeny will include only branches above the first split.

  \code{bisse.tree} may return an extinct phylogeny, and
  \code{bisse.trees} might return extinct phylogenies if
  \code{include.extinct} is \code{TRUE}.
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

\author{Richard G. FitzJohn}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(3)
phy <- tree.bisse(pars, max.taxa=30, x0=0)
phy$tip.state

plot(phy, label.offset=1, cex=.8)
tiplabels(pch=21, bg=phy$tip.state+1, adj=1)
nodelabels(pch=21, bg=phy$node.state+1)
}

\keyword{models}
