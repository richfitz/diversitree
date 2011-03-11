\name{make.mkn.deriv}
\alias{make.mkn.deriv}

\title{Mk2 and Mk-n Models of character evolution [temporary extension
  with derivatives]}

\description{Currently a restricted-feature-set version of
  \code{\link{make.mkn}}, but including gradient (first derivative)
  information along with the likelihood.  Currently lacks several
  features of the full \code{mkn} model: no support for missing trait
  data, no support for different root treatments, no optimised 2-state
  method.

  In the future, this will roll back into \code{\link{make.mkn}}, so do
  not use or depend on this!
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be
  integers from 1 to \code{k}.}

\item{k}{Number of states to model.}
}

\seealso{
  \code{\link{make.mkn}}, which is more general, and from which this
  inherits support methods for ML and MCMC inference.
}

\author{Richard G. FitzJohn}
\keyword{models}





