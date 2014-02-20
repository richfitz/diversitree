\name{make.geosse.split}
\alias{make.geosse.split}
\alias{make.geosse.uneven}

\title{Geographic State Speciation and Extinction Model: Split Models}

\description{Create a likelihood function for a GeoSSE model where the
  tree is partitioned into regions with different parameters.}

\usage{
make.geosse.split(tree, states, nodes, split.t,
                  sampling.f=NULL, strict=TRUE, control=list())
make.geosse.uneven(tree, states, nodes, split.t,
                   sampling.f=NULL, strict=TRUE, control=list())
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be an
    integer between 0 and 2: see \code{\link{make.geosse}}.  This vector
    must have names that correspond to the tip labels in the
    phylogenetic tree (\code{tree$tip.label}).  Unresolved clades are
    not supported.}
  
  \item{nodes}{Vector of nodes that will be split (see Details).}
  
  \item{split.t}{Vector of split times, same length as \code{nodes} (see
    Details).}
  
  \item{sampling.f}{Vector of length 3 where \code{sampling.f[i]}
    is the proportion of species in state \code{i} that are present in
    the phylogeny.  A value of \code{c(0.5, 0.75, 1)} means that half of
    species in state 0, three quarters of species in state 1, and all
    species in state 2 are included in the phylogeny.  By default all
    species are assumed to be known.}

  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are integers on \code{0:2}.  If \code{strict} is
    \code{TRUE} (the default), then the additional check is made that
    \emph{every} state is present.  The likelihood models tend to be
    poorly behaved where states are missing, but there are cases
    (missing intermediate states for meristic characters) where allowing
    such models may be useful.}

  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}
}

\details{Branching times can be controlled with the \code{split.t}
argument.  If this is \code{Inf}, split at the base of the branch (as in
MEDUSA).  If \code{0}, split at the top (closest to the present, as in
the new option for MEDUSA).  If \code{0 < split.t < Inf} then we split
at that time on the tree (zero is the present, with time growing
backwards).

The \code{nodes} at the top of the split location are specified as a
vector of node names.  For example, a value of \code{c("nd10",
"nd12")} means that the splits are along the branches leading from each
of these nodes towards the root.
}

\author{Emma E. Goldberg}
\keyword{models}
