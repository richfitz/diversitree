\name{utilities}
\alias{get.descendants}
\title{Utility Functions}
\description{
  These are utility functions that are used internally by diversitree,
  but which might be more generally useful.

  Currently only \code{get.descendants} docuemnted here, which
  determines which species or nodes are descended from a particular
  node.

}

\usage{
get.descendants(node, tree, tips.only=FALSE, edge.index=FALSE)
}

\arguments{
  \item{node}{A node, either a name in \code{tree$node.label}, an
    integer in \code{1..tree$Nnode}, or in
    \code{length(tree$tip.label)..(length(tree$tip.label)+tree$Nnode)}.}

  \item{tree}{A phylogenetic tree, in ape's \code{phylo} format.}

  \item{tips.only}{Logical: return only descendant indices of tip
    species?}

  \item{edge.index}{Logical: return the row indices in the edge matrix?}
}

\author{Richard G. FitzJohn}

\keyword{utilities}


