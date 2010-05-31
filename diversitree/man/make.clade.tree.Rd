\name{make.clade.tree}
\alias{make.clade.tree}
\alias{clades.from.polytomies}
\alias{clades.from.classification}

\title{Make a "Clade Tree"}

\description{This function makes a \dQuote{clade tree}, where tips
  represent clades.  It is designed to make working with unresolved
  clade information in \code{\link{make.bisse}} more straightforward.
  \code{clade.tree} objects have their own plotting methods.}

\usage{
make.clade.tree(tree, clades)
clades.from.polytomies(tree)
clades.from.classification(tree, class, check=TRUE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{phylo} (Ape) format.}
  \item{clades}{A list, where the name of each element represents a tip
    in \code{tree} and each element is a character vector containing the
    names of species contained within that clade.}
  \item{class}{A vector along \code{tree$tip.label}, indicating the
    class to which each tip belongs.}
  \item{check}{Logical, indicating whether a (potentially slow) check
    will be done to ensure that all resulting clades are reciprocally
    monophyletic within the tree.}
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
\author{Richard G. FitzJohn}
