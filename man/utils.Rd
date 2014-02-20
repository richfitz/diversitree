\name{utilities}
\alias{get.descendants}
\alias{run.cached}
\alias{expand.parameters}
\alias{drop.likelihood}
\alias{get.likelihood}
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
run.cached(filename, expr, regenerate=FALSE)
expand.parameters(p, lik.new, repl=0, target=argnames(lik.new))
get.likelihood(object)
drop.likelihood(object)
}

\arguments{
  \item{node}{A node, either a name in \code{tree$node.label}, an
    integer in \code{1..tree$Nnode}, or in
    \code{length(tree$tip.label)..(length(tree$tip.label)+tree$Nnode)}.}

  \item{tree}{A phylogenetic tree, in ape's \code{phylo} format.}

  \item{tips.only}{Logical: return only descendant indices of tip
    species?}

  \item{edge.index}{Logical: return the row indices in the edge matrix?}

  \item{filename}{Name of the file to store cached results}

  \item{expr}{Expression to evaluate}

  \item{regenerate}{Logical: force re-evaluation of expr and
    regeneration of filename?}

  \item{object}{For \code{drop.likelihood}, an object that has a
    \code{likelihood} attribute to be removed (saves space on object
    save); for \code{get.likelihood}, retrieves the function.}

  \item{p, lik.new, repl, target}{Undocumented currently}
}

\author{Richard G. FitzJohn}

\keyword{utilities}
