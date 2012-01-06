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
make.bm(tree, states, states.sd=0, method="vcv")
make.ou(tree, states, states.sd=0)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of continuous valued character states.  This
    vector must be named with the tip labels of \code{tree}.}
  
  \item{states.sd}{An optional vector of measurement errors, as standard
    deviation around the mean.  If a single value is given it is used
    for all tips, otherwise the vector must be named as for
    \code{states}.}

  \item{method}{A character string indicating if the likelihood should
    be computed using the VCV matrix (\code{method="vcv"}, the default,
    and the approach used by Geiger), or directly using along-branch
    calculations (\code{method="direct"}, more similar to how most
    other models are computed with diversitree).}
}

\details{
  This function is heavily based on \code{fitContinuous} in the
  \code{geiger} package.  However, this version only fits Brownian
  motion models.
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

## With the direct calculations
fit2 <- find.mle(make.bm(bird.orders, x, method="direct"), .1)

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
