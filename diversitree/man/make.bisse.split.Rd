\name{make.bisse.split}
\alias{make.bisse.split}

\title{Binary State Speciation and Extinction Model: Split Models}

\description{Create a likelihood function for a BiSSE model where the
  tree is partitioned into regions with different parameters.}

\usage{
make.bisse.split(tree, states, nodes, split.t, unresolved=NULL,
                 sampling.f=NULL, nt.extra=10, strict=TRUE, control=list())
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be 0 or
    1, or \code{NA} if the state is unknown.  This vector must have
    names that correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips
    corresponding to unresolved clades, the state should be \code{NA}.}
  
  \item{nodes}{Vector of nodes that will be split (see Details).}
  
  \item{split.t}{Vector of split times, same length as \code{nodes} (see
    Details).}
  
  \item{unresolved}{Unresolved clade information: see section below for
    structure.}
  
  \item{sampling.f}{Vector of length 2 with the estimated proportion of
    extant species in state 0 and 1 that are included in the phylogeny.
    A value of \code{c(0.5, 0.75)} means that half of species in state 0
    and three quarters of species in state 1 are included in the
    phylogeny.  By default all species are assumed to be known.}
  
  \item{nt.extra}{The number of species modelled in unresolved clades
    (this is in addition to the largest observed clade).}

  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are 0 and 1 only.  If \code{strict} is \code{TRUE}
    (the default), then the additional check is made that \emph{every}
    state is present.  The likelihood models tend to be poorly behaved
    where states are missing.}
  
  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}
}


\details{Branching times can be controlled with the \code{split.t}
argument.  If this is \code{Inf}, split at the base of the branch (as in
MEDUSA).  If \code{0}, split at the top (closest to the present, as in
the new option for MEDUSA).  If \code{0 < split.t < Inf} then we split
at that time on the tree (zero is the present, with time growing
backwards).

TODO: Describe \code{nodes} and \code{split.t} here.
}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(546)
phy <- tree.bisse(pars, max.taxa=30, x0=0)

## Here is the phylogeny:
plot(phy, show.node.label=TRUE, label.offset=.1, font=1, cex=.75,
     no.margin=TRUE)

## Here is a plain BiSSE function for comparison:
lik.b <- make.bisse(phy, phy$tip.state)
lik.b(pars) # -93.62479

## Split this phylogeny at three points: nd15, nd18 and nd26
nodes <- c("nd15", "nd18", "nd26")

## This is the index in ape's node indexing:
nodes.i <- match(nodes, phy$node.label) + length(phy$tip.label)

nodelabels(node=nodes.i, pch=19, cex=2, col="#FF000099")

## Split these at the base of the branch
split.t <- Inf

## To make a split BiSSE function, pass the node locations and times in:
lik.s <- make.bisse.split(phy, phy$tip.state, nodes.i, split.t)

## The parameters must be a list of the same length as the number of
## partitions.  Partition '1' is the root partition, and partition i is
## the partition rooted at the node[i-1]
pars4 <- rep(list(pars), 4)
pars4

## Run the likelihod calculation:
lik.s(pars4) # -93.62479

## These are basically identical (to acceptable tolerance)
lik.s(pars4) - lik.b(pars)

## You can use the labelled nodes rather than indices:
lik.s2 <- make.bisse.split(phy, phy$tip.state, nodes, split.t)
identical(lik.s(pars4), lik.s2(pars4))

## This also works where some tips are unresolved clades.  Here are a
## few:
unresolved <-
  data.frame(tip.label=c("sp12", "sp32", "sp9", "sp22", "sp11"),
             Nc=c(2,5,3,2,5), n0=c(1, 4, 3, 2, 4), n1=c(1, 1, 0, 0, 1))

## Plain BiSSE with unresolved clades:
lik.u.b <- make.bisse(phy, phy$tip.state, unresolved=unresolved)
lik.u.b(pars) # -139.3669

## Split BiSSE with unresolved clades:
lik.u.s <- make.bisse.split(phy, phy$tip.state, nodes, split.t,
                            unresolved=unresolved)
lik.u.s(pars4) # -139.3669

lik.u.b(pars) - lik.u.s(pars4) # numerical error only
}

\author{Richard G. FitzJohn}
\keyword{models}
