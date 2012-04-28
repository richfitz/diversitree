\name{make.quasse.split}
\alias{make.quasse.split}

\title{Quantitative State Speciation and Extinction Model: Split Models}

\description{Create a likelihood function for a QuaSSE model where the
  tree is partitioned into regions with different parameters.}

\usage{
make.quasse.split(tree, states, states.sd, lambda, mu, nodes, split.t,
                  control=NULL, sampling.f=NULL)
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be a
    numeric real values.  Missing values (\code{NA}) are not yet handled.
    This vector must have names that correspond to the tip labels in the
    phylogenetic tree  (\code{tree$tip.label}).}

  \item{states.sd}{A scalar or vector corresponding to the standard
    error around the mean in states (the initial probability
    distribution is assumed to be normal).}

  \item{lambda}{A function to use as the speciation function.  The first
    argument of this must be \code{x} (see Details).}

  \item{mu}{A function to use as the extinction function.  The first
    argument of this must be \code{x} (see Details.)}

  \item{nodes}{Vector of nodes that will be split (see Details).}
  
  \item{split.t}{Vector of split times, same length as \code{nodes} (see
    Details).}

  \item{control}{A list of parameters for tuning the performance of the
    integrator.  A guess at reasonble values will be made here.  See
    Details in \code{\link{make.quasse}} for possible entries.}

  \item{sampling.f}{Scalar with the estimated proportion of extant
    species that are included in the phylogeny.  A value of \code{0.75}
    means that three quarters of extant species are included in the
    phylogeny.  By default all species are assumed to be known.}
}

\details{Branching times can be controlled with the \code{split.t}
argument.  If this is \code{Inf}, split at the base of the branch (as in
MEDUSA).  If \code{0}, split at the top (closest to the present, as in
the new option for MEDUSA).  If \code{0 < split.t < Inf} then we split
at that time on the tree (zero is the present, with time growing
backwards).

TODO: Describe \code{nodes} and \code{split.t} here.
}

\author{Richard G. FitzJohn}
\keyword{models}
