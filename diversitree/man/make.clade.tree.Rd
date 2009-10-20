\name{make.clade.tree}
\alias{make.clade.tree}
\alias{polytomies.to.clades}

\title{Make a "Clade Tree"}

\description{This function makes a \dQuote{clade tree}, where tips
  represent clades.  It is designed to make working with unresolved
  clade information in \code{\link{make.bisse}} more straightforward.
  Soon, \code{clade.tree} objects will have their own plotting methods,
  though.}

\usage{
make.clade.tree(tree, clades)
polytomies.to.clades
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{phylo} (Ape) format.}
  \item{clades}{A list, where the name of each element represents a tip
    in \code{tree} and each element is a character vector containing the
    names of species contained within that clade.}
}

\details{
  The idea here is that \code{make.bisse} takes a tree and a named
  character state vector.  If the phylogenetic tree contains information
  about the membership of clades, then the unresolved clade information
  can be constructed automatically.  The names chosen should therefore
  reflect the names used in the state information.
}

\note{
  Currently, \code{clade.tree} objects work poorly with some \code{ape}
  functions.
}

\keyword{manip}
