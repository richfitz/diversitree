\name{make.bm}
\alias{make.bm}
\alias{make.ou}

\title{Brownian Motion and Ornstein-Uhlenbeck Character Evolution}

\description{Create a likelihood function for models of simple Brownian
  Motion or Ornstein-Uhlenbeck (OU) character evolution.  This function
  creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.bm(tree, states, states.sd=0, control=list())
make.ou(tree, states, states.sd=0, control=list())
}

\arguments{
  \item{tree}{A bifurcating phylogenetic tree, in \code{ape}
    \dQuote{phylo} format.}
  
  \item{states}{A vector of continuous valued character states.  This
    vector must be named with the tip labels of \code{tree}.}
  
  \item{states.sd}{An optional vector of measurement errors, as standard
    deviation around the mean.  If a single value is given it is used
    for all tips, otherwise the vector must be named as for
    \code{states}.}

  \item{control}{A list of control parameters.  See details below.}
}

\details{The \code{control} argument is a named list of options.

  The main option is \code{method}.  Specifying
  \code{control=list(method="vcv")} uses a variance-covariance matrix
  based approach to compute the likelihood.  This is similar to the
  approach used by geiger, and is the default for \code{make.bm} (but is
  not available for \code{make.ou}).  The alternative is to specify
  \code{control=list(method="pruning")} which uses the transition density
  function for brownian motion along each branch, similar to how most
  methods in diversitree are computed.  This second approach is much
  faster for very large trees.

  When \code{method="vcv"} is specified, \code{backend="R"} or
  \code{backend="C"} may also be provided, which switch between a slow
  (and stable) R calculator and a fast (but less extensively tested) C
  calculator.  \code{backend="R"} is currently the default.

  The VCV-based \code{make.bm} function is heavily based on
  \code{fitContinuous} in the \code{geiger} package.

  For non-ultrametric trees with OU models, computed likelihoods will
  differ because of the different root treatments.

  % TODO: Describe root treatments!
}

\seealso{
  See
  \url{http://www.zoology.ubc.ca/prog/diversitree/examples/ou-nonultrametric/}
  for a discussion about calculations on non-ultrametric trees.
}

% TODO: Add example section and show equivalence to geiger?
\examples{
## Random data (following APE)
data(bird.orders)
set.seed(1)
x <- structure(rnorm(length(bird.orders$tip.label)),
               names=bird.orders$tip.label)

## With the VCV approach
fit1 <- find.mle(make.bm(bird.orders, x), .1)

## With the pruning calculations
lik.pruning <- make.bm(bird.orders, x, control=list(method="pruning"))
fit2 <- find.mle(lik.pruning, .1)

## All the same (need to drop the function from this though)
all.equal(fit1[-7], fit2[-7])

## If this is the same as the estimates from Geiger, to within the
## tolerances expected for the calculation and optimisation:
\dontrun{
fit3 <- fitContinuous(bird.orders, x)
all.equal(fit3$Trait1$lnl, fit1$lnLik)
all.equal(fit3$Trait1$beta, fit1$par, check.attributes=FALSE)
}

}

\author{Richard G. FitzJohn}
\keyword{models}
