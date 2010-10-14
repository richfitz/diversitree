\name{make.bd.split}
\alias{make.bd.split}

\title{Constant Rate Birth-Death Models: Split Models}

\description{Create a likelihood function for a birth-death model where
  the tree is partitioned into regions with different parameters.}

\usage{
make.bd.split(tree, nodes, split.t, sampling.f=NULL, unresolved=NULL)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{nodes}{Vector of nodes that will be split (see Details).}
  \item{split.t}{Vector of split times, same length as \code{nodes} (see
    Details).}
  \item{sampling.f}{Probability of an extant species being included in
    the phylogeny (sampling fraction).  By default, all extant species are
    assumed to be included.}
  \item{unresolved}{Unresolved clade information.  This is a named
    vector, with the number of species as the value and names
    corresponding to tip labels.  Tips that represent a single species
    should not be included in this vector.  For example
    \code{sp1=10, sp2=2}, would mean that \code{sp1} represents 10
    species, while \code{sp2} represents two.  These labels must exist
    in \code{tree$tip.label} and all other tips are assumed to
    represent one species.}
}

\details{Branching times can be controlled with the \code{split.t}
argument.  If this is \code{Inf}, split at the base of the branch (as in
MEDUSA).  If \code{0}, split at the top (closest to the present, as in
the new option for MEDUSA).  If \code{0 < split.t < Inf} then we split
at that time on the tree (zero is the present, with time growing
backwards).

TODO: Describe \code{nodes} and \code{split.t} here.

This function is related to MEDUSA (Alfaro et al. 2009), but does not
include any of the code for efficiently moving between different splits
(split creation here is fairly slow).  The primary use for this model is
for generating starting points for state dependent split models (e.g.,
\code{\link{make.bisse.split}}) and testing \emph{a priori} splits.
}

\examples{
set.seed(1)
pars <- c(.1, .03)
phy <- trees(pars, "bd", max.taxa=30)[[1]]

## Here is the phylogeny:
plot(phy, show.node.label=TRUE, label.offset=.1, font=1, cex=.75,
     no.margin=TRUE)

## Construct the plain likelihood function as a benchmark:
lik <- make.bd(phy)
lik(pars) # -21.74554

## Split this phylogeny at three points: nd11, nd13 and nd26
nodes <- c("nd11", "nd13", "nd26")

## This is the index in ape's node indexing:
nodes.i <- match(nodes, phy$node.label) + length(phy$tip.label)

nodelabels(node=nodes.i, pch=19, cex=2, col="#FF000099")

## To make a split likelihood function, pass the node locations and times in:
lik.s <- make.bd.split(phy, nodes.i)

## The parameters must be a list of the same length as the number of
## partitions.  Partition '1' is the root partition, and partition i is
## the partition rooted at the node[i-1]
pars4 <- rep(pars, 4)
names(pars4) <- argnames(lik.s)

## Run the likelihod calculation:
lik.s(pars4) # -21.74554

## These are basically identical (to acceptable tolerance)
lik.s(pars4) - lik(pars)

## You can use the labelled nodes rather than indices:
lik.s2 <- make.bd.split(phy, nodes)
identical(lik.s(pars4), lik.s2(pars4))

## All the usual ML/MCMC functions work as before:
fit <- find.mle(lik.s, pars4)
}

\author{Richard G. FitzJohn}
\keyword{models}
