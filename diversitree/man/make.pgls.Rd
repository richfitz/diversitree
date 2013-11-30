\name{make.pgls}
\alias{make.pgls}

\title{Phylogenetic Generalised Least Squares}

\description{Generate the likelihood function that underlies PGLS
  (Phylogenetic Generalised Least Squares).  This is a bit of a misnomer
  here, as you may not be interested in least squares (e.g., if using
  this with \code{\link{mcmc}} for Bayesian inference).}

\usage{
make.pgls(tree, formula, data, control=list())
}

\arguments{
  \item{tree}{A bifurcating phylogenetic tree, in \code{ape}
    \dQuote{phylo} format.}

  \item{formula}{A model formula; see \code{\link{lm}} for details on
    formulae; the interface is the same here.}

  \item{data}{A data frame containing the variables in the model.  If
    not found in \code{data}, the variables are taken from
    \code{environment(formula)}, typically the environment from
    which this function is called.  That may perform badly with
    reconciling with species names, however.}

  \item{control}{A list of control parameters.  Currently the only
    option is the key \dQuote{method} which can be \code{"vcv"} for the
    traditional variance-covariance approach (slow for large trees) or
    \code{"contrasts"} for the contrasts-based approach outlined in
    Freckleton (2012).}
}

\references{
  Freckleton R.P. 2012. Fast likelihood calculations for comparative
  analyses. Methods in Ecology and Evolution 3: 940-947.
}

\author{Richard G. FitzJohn}
\keyword{models}
