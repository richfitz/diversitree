\name{make.bm}
\alias{make.bm}

\title{Brownian Motion Character Evolution}

\description{Create a likelihood function for models of simple Brownian
  Motion character evolution.  This function creates a likelihood
  function that can be used in \link[=find.mle]{maximum likelihood} or
  \link[=mcmc]{Bayesian} inference.}

\usage{
make.bm(tree, states, meserr=NULL)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of continuous valued character states.  This
    vector must be named with the tip labels of \code{tree}.}
  
  \item{meserr}{An optional vector of measurement errors, as standard
    deviation around the mean.  If a single value is given it is used
    for all tips, otherwise the vector must be named as for
    \code{states}.}
}

\details{
  This function is heavily based on \code{fitContinuous} in the
  \code{geiger} package.  However, this version only fits Brownian
  motion models.
}

% TODO: Add example section and show equivalence to geiger?
%\examples{
%}

\author{Richard G. FitzJohn}
\keyword{models}



